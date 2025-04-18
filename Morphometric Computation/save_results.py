import pandas as pd
import os

def save_results(mouse_id: str, sex: str, genotype: str, age: str, 
                pixel_to_micrometer, physical_width, physical_height,
                side_stats: dict, sub_geometry: dict):
    """
    results are saved in 34 CSV files
    1. summary_stats.csv - include sides (proximal, distal) and global stat
    2. distribution_data.csv - the area distribution for each plaque and each neuron
    3. geometry_metadata.csv - data about the subiculum
    4. spatial_data.csv - spatial distribution of centroids of each plaque and neuron
    """
    
    # ================== 1. summary_stats ==================
    summary_data = []
    for side in ['left', 'right', 'total']:
        # SIDE
        stats = side_stats.get(side, {})
        sub_area = sub_geometry['geometry']['area_um2']
        
        # normalize total area using subiculum area
        cell_area_pct = (stats.get('cell_area', 0) / sub_area) if sub_area > 0 else 0
        plaque_area_pct = (stats.get('plaque_area', 0) / sub_area) if sub_area > 0 else 0
        intra_pct = (stats.get('intra_plaque', 0) / sub_area) if sub_area > 0 else 0
        
        row = {
            'mouse_id': mouse_id,
            'sex': sex,
            'genotype': genotype,
            'age': age,
            'pixel_to_micrometer': pixel_to_micrometer,
            'width_in_micrometer': physical_width,
            'height_in_micrometer': physical_height,
            'side': side.upper(),
            
            # CELL
            'n_cells': stats.get('cell_count', 0),
            'total_cell_area': stats.get('cell_area', 0),
            'mean_cell_area': stats.get('cell_area', 0)/stats.get('cell_count', 1) if stats.get('cell_count', 0) > 0 else 0,
            'cell_area_pct': cell_area_pct,
            
            # AB
            'n_plaques': stats.get('plaques', 0),
            'total_plaque_area': stats.get('plaque_area', 0),
            'mean_plaque_area': stats.get('plaque_area', 0)/stats.get('plaques', 1) if stats.get('plaques', 0) > 0 else 0,
            'plaque_area_pct': plaque_area_pct,
            'intra_plaque_pct': intra_pct,
            
            # GLOBAL
            'sub_area': sub_area,
            'sub_width': sub_geometry['geometry']['width_um'],
            'sub_height': sub_geometry['geometry']['height_um']
        }
        summary_data.append(row)
    
    _save_to_csv(pd.DataFrame(summary_data), 'summary_stats.csv')
    
    # ================== 2. Save Distribution ==================
    dist_data = []
    
    # global cell area distribution
    for area in side_stats.get('cell_global_dist', []):
        dist_data.append({
            'mouse_id': mouse_id,
            'type': 'cell',
            'side': 'GLOBAL',
            'area_um2': area
        })
    # global ab area distrbution
    for area in side_stats.get('ab_global_dist', []):
        dist_data.append({
            'mouse_id': mouse_id,
            'type': 'plaque',
            'side': 'GLOBAL',
            'area_um2': area
        })
    
    # sides
    for side in ['left', 'right']:
        for area in side_stats.get(side, {}).get('cell_dist', []):
            dist_data.append({
                'mouse_id': mouse_id,
                'type': 'cell',
                'side': side.upper(),
                'area_um2': area
            })
            
        for area in side_stats.get(side, {}).get('ab_dist', []):
            dist_data.append({
                'mouse_id': mouse_id,
                'type': 'plaque',
                'side': side.upper(),
                'area_um2': area
            })
    
    if dist_data:
        _save_to_csv(pd.DataFrame(dist_data), 'distribution_data.csv')
    
    # ================== 3. subiculum geometry ==================
    geo_data = {
        'mouse_id': mouse_id,
        'sub_area_um2': sub_geometry['geometry']['area_um2'],
        'sub_width_um': sub_geometry['geometry']['width_um'],
        'sub_height_um': sub_geometry['geometry']['height_um'],
        
        # midline connecting two ends
        'midline_slope': sub_geometry['midline'].get('slope'),
        'midline_intercept': sub_geometry['midline'].get('intercept'),
        'midline_midpoint_x': sub_geometry['midline'].get('midpoint_um', (None, None))[0],
        'midline_midpoint_y': sub_geometry['midline'].get('midpoint_um', (None, None))[1],
        'midline_is_vertical': sub_geometry['midline'].get('is_vertical', False),
        'midline_is_horizontal': sub_geometry['midline'].get('is_horizontal', False),
        
        # dividing line, perpendicular to midline
        'divline_slope': sub_geometry['dividing_line'].get('slope'),
        'divline_intercept': sub_geometry['dividing_line'].get('intercept'),
        'divline_is_vertical': sub_geometry['dividing_line'].get('is_vertical', False),
        'divline_is_horizontal': sub_geometry['dividing_line'].get('is_horizontal', False)
    }
    
    _save_to_csv(pd.DataFrame([geo_data]), 'geometry_metadata.csv')

    # ================== 4. Centroids distribution ==================
    spatial_data = []
    
    # Cell centroids
    for side in ['left', 'right']:
        stats = side_stats.get(side, {})
        for idx, (x, y) in enumerate(stats.get('cell_centroids', [])):
            spatial_data.append({
                'mouse_id': mouse_id,
                'type': 'cell',
                'side': side.upper(),
                'id': idx + 1,  # Cell ID
                'x_um': x,
                'y_um': y
            })
    
    # AB centroids
    for side in ['left', 'right']:
        stats = side_stats.get(side, {})
        for idx, (x, y) in enumerate(stats.get('ab_centroids', [])):
            spatial_data.append({
                'mouse_id': mouse_id,
                'type': 'plaque',
                'side': side.upper(),
                'id': idx + 1,  # AB ID
                'x_um': x,
                'y_um': y
            })
    
    if spatial_data:
        _save_to_csv(pd.DataFrame(spatial_data), 'spatial_data.csv')

def _save_to_csv(df: pd.DataFrame, filename: str):
    file_exists = os.path.exists(filename)
    df.to_csv(filename, mode='a', index=False, header=not file_exists)