# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 10:23:03 2025

@author: mandrolk1
"""


import numpy as np
from scipy.spatial import cKDTree
import networkx as nx

def parse_first_timestep(input_filename):
    """Парсить перший таймстеп із файлу, повертає лише атомні дані"""
    atoms = []
    with open(input_filename, 'r') as infile:
        inside_timestep = False
        in_atoms_section = False
        for line in infile:
            line = line.strip()
            if line.startswith("ITEM: TIMESTEP"):
                if not inside_timestep:
                    inside_timestep = True
                else:
                    break  # припиняємо після першого таймстепу
                continue
            if inside_timestep:
                if line.startswith("ITEM: ATOMS"):
                    in_atoms_section = True
                    continue
                if in_atoms_section and line:
                    parts = line.split()
                    if len(parts) >= 6:  # перевірка наявності id, type, name, x, y, z
                        atoms.append(parts[:6])  # зберігаємо лише потрібні стовпці
    return np.array(atoms, dtype=str) if atoms else np.empty((0, 6))

def compute_distance(o_coord, h_coord, box_size):
    """
    Обчислює 3D відстань між оксигеном та гідрогеном з урахуванням
    періодичних умов у x та y.
    """
    delta = o_coord - h_coord  # різниця по x,y,z
    # Коригуємо x та y з урахуванням PBC
    delta[0] = delta[0] - box_size[0] * np.round(delta[0] / box_size[0])
    delta[1] = delta[1] - box_size[1] * np.round(delta[1] / box_size[1])
    return np.linalg.norm(delta)

def get_water_molecules(input_filename, output_file, bond_cutoff=1.2):
    """Ідентифікує молекули води за допомогою побудови графу та максимального парного зіставлення"""
    # Розміри коробки: періодичність застосовується у x та y, по z – немає
    box_size = np.array([60.34, 139.34, np.inf])
    
    # Парсимо дані першого таймстепу
    atom_data = parse_first_timestep(input_filename)
    if atom_data.size == 0:
        print("Дані не знайдено в першому таймстепі!")
        return
    
    # Отримуємо id, назви та координати
    atom_ids = atom_data[:, 0].astype(int)
    atom_names = atom_data[:, 2]
    coords = atom_data[:, 3:6].astype(float)
    
    # Розділяємо на оксигени (O) та гідрогени (H)
    o_mask = (atom_names == 'O')
    h_mask = (atom_names == 'H')
    
    o_ids = atom_ids[o_mask]
    o_coords = coords[o_mask]
    h_ids = atom_ids[h_mask]
    h_coords = coords[h_mask]
    
    if len(o_ids) == 0 or len(h_ids) < 2:
        print("Недостатньо атомів для формування молекул води!")
        return
    
    # Створюємо дерево для гідрогенів використовуючи лише x та y координати
    h_coords_xy = h_coords[:, :2]
    box_xy = box_size[:2].tolist()
    h_tree = cKDTree(h_coords_xy, boxsize=box_xy)
    
    # Створимо граф, де кожен оксиген буде представленим двома вершинами,
    # адже він повинен мати два зв'язки. Для унікальності ідентифікації використовуємо tuple:
    # (o_id, 'a') та (o_id, 'b').
    G = nx.Graph()
    oxygen_nodes = []
    for idx, o_id in enumerate(o_ids):
        node_a = (o_id, 'a')
        node_b = (o_id, 'b')
        oxygen_nodes.extend([node_a, node_b])
        # Додаємо вершини для оксигену до графу
        G.add_node(node_a, bipartite=0, o_index=idx, o_id=o_id)
        G.add_node(node_b, bipartite=0, o_index=idx, o_id=o_id)
    
    # Додаємо вершини для гідрогенів
    hydrogen_nodes = []
    for j, h_id in enumerate(h_ids):
        h_node = ("H", int(h_id))  # унікальний ідентифікатор для гідрогена
        hydrogen_nodes.append(h_node)
        G.add_node(h_node, bipartite=1, h_index=j, h_id=h_id)
    
    # Для кожного оксигену знаходимо кандидатів-сусідів серед гідрогенів
    # Скористаємося пошуком по 2D (x,y) для попереднього відбору, а потім перевіримо повну 3D відстань.
    for idx, (o_id, o_coord) in enumerate(zip(o_ids, o_coords)):
        # Отримуємо x,y оксигену
        o_xy = o_coord[:2]
        # Пошук у дереві – query_ball_point повертає список індексів кандидатів
        candidate_indices = h_tree.query_ball_point(o_xy, r=bond_cutoff)
        
        # Для кожного кандидата перевіряємо повну відстань з урахуванням PBC
        for j in candidate_indices:
            h_coord = h_coords[j]
            r = compute_distance(o_coord, h_coord, box_size)
            if r < bond_cutoff:
                # Додаємо ребро від обох копій оксигену до цього гідрогену
                h_node = ("H", int(h_ids[j]))
                node_a = (o_id, 'a')
                node_b = (o_id, 'b')
                # Додаємо ребро, якщо його ще нема (можна задати вагу або відстань, якщо потрібно оптимізувати)
                G.add_edge(node_a, h_node, weight=r)
                G.add_edge(node_b, h_node, weight=r)
    
    # Виконуємо максимальне парне зіставлення (двостороннє – Hopcroft-Karp)
    matching = nx.algorithms.bipartite.matching.hopcroft_karp_matching(G, top_nodes=oxygen_nodes)
    
    # Збираємо результати: для кожного оксигену, якщо обидві копії
    # успішно зіставлені з гідрогенами, вважаємо, що знайдена молекула води.
    water_molecules = []
    # Оскільки matching містить пари з обох сторін, обробляємо лише вершини оксигену (oxygen_nodes)
    processed_oxygens = set()
    for node in oxygen_nodes:
        if node in matching:
            o_id = node[0]
            # Щоб не дублювати, перевіряємо, чи вже обробили цей оксиген
            if o_id in processed_oxygens:
                continue
            # Отримуємо обидві копії для даного оксигену
            node_a = (o_id, 'a')
            node_b = (o_id, 'b')
            # Якщо обидва присутні в зіставленні, отримуємо їхніх сусідів
            if node_a in matching and node_b in matching:
                h1 = matching[node_a]
                h2 = matching[node_b]
                # Переконуємося, що отримані вершини дійсно є гідрогенами
                if isinstance(h1, tuple) and h1[0] == "H" and isinstance(h2, tuple) and h2[0] == "H":
                    water_molecules.append((o_id, h1[1], h2[1]))
                    processed_oxygens.add(o_id)
    
    # Зберігаємо результати у файл
    with open(output_file, 'w') as f:
        f.write("O_id\tH1_id\tH2_id\n")
        for mol in water_molecules:
            f.write(f"{mol[0]}\t{mol[1]}\t{mol[2]}\n")
    
    #print(f"Знайдено {len(water_molecules)} молекул води. Результати збережено у {output_file}")

# Приклад використання
input_file = "dynamics.xyz"
output_file = "water_molecules_ids.dat"
bond_cutoff = 1.2

get_water_molecules(input_file, output_file, bond_cutoff)
