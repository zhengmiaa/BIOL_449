import pandas as pd
import os

# Save the results to a CSV file
def save_results(cell_mask_area, unique_cell_masks, cell_mask_area_dist, ab_mask_area, unique_ab_masks, ab_mask_area_dist, intra_overlap_area_um2, 
                 sub_area, sub_width_um, sub_height_um, sub_center_um, 
                 cell_centers, ab_centers, min_distances, 
                 age, pixel_to_micrometer, physical_width, physical_height):
    mouse_id = input("Enter Mouse ID: ")

    sex_choice = input("Choose Sex (1. Female, 2.Male): ")
    sex_mapping = {'1': 'F', '2':'M'}
    sex = sex_mapping.get(sex_choice, "NA")

    genotype_choice = input("Choose Genotype (1.WT, 2.5XFAD): ")
    genotype_mapping = {'1': 'WT', '2':'5XFAD'}
    genotype = genotype_mapping.get(genotype_choice, "NA")
    
    # age_choice = input("Choose Age (1. 6 weeks, 2. 10 weeks, 3. 6 months):")
    # age_mapping = {'1': '6 weeks','2': '10 weeks','3': '6 months'}
    # age = age_mapping.get(age_choice, "NA")

    average_cell_size = cell_mask_area / unique_cell_masks if unique_cell_masks > 0 else 0  # Avoid division by zero
    percent_cell_area = cell_mask_area / sub_area if sub_area > 0 else 0

    average_ab_size = ab_mask_area / unique_ab_masks if unique_ab_masks > 0 else 0  # Avoid division by zero
    percent_ab_area = ab_mask_area / sub_area if sub_area > 0 else 0

    percent_ab_area_intra = intra_overlap_area_um2 / sub_area if sub_area > 0 else 0
    data = {
            'Mouse_ID': mouse_id,
            'Sex': sex,
            'Genotype': genotype,
            'Age': age,
            
            'um/Pixel_Factor': pixel_to_micrometer,
            'Pic_Width': physical_width,
            'Pic_Height': physical_height,
            
            'SUB_Area_um2': sub_area,
            'SUB_Width_um': sub_width_um, 
            'SUB_Height_um': sub_height_um, 
            'SUB_Center_Coord_um': sub_center_um, 
            'Cell_Count': unique_cell_masks,

            'Avg_Cell_Size_um2': average_cell_size,
            'Percent_Cell_Area_in_SUB': percent_cell_area,
            'Cell_Size_Dist__um2': cell_mask_area_dist,
            
            'AB_Count': unique_ab_masks,
            'Avg_AB_Size_um2': average_ab_size,
            'Percent_AB_Area_in_SUB': percent_ab_area,
            'Percent_Intracellular_AB_Area_in_SUB': percent_ab_area_intra,
            'AB_Size_Dist_um2': ab_mask_area_dist,

            'Cell_Centroids_um': cell_centers,
            'AB_Centroids_um': ab_centers,
            'Min_Distances_Cell_to_AB_um)': min_distances
        }
    
    # Create a DataFrame and save to CSV
    csv_file_path = 'segmentation_results_all.csv'
    file_exists = os.path.exists(csv_file_path)
    df = pd.DataFrame([data])
    df.to_csv(csv_file_path, mode='a', index=False, header=not file_exists)  # Append if exists
    print(f'Results saved.')



import pandas as pd
import os

def save_results(mouse_id: str, sex: str, genotype: str, age: str, 
                pixel_to_micrometer, physical_width, physical_height,
                side_stats: dict, sub_geometry: dict):
    """
    保存结果为三个CSV文件：
    1. summary_stats.csv - 分侧统计摘要
    2. distribution_data.csv - 分侧面积分布
    3. geometry_metadata.csv - 几何参数元数据
    """
    
    # ================== 1. 保存统计摘要 ==================
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
    """统一保存函数"""
    file_exists = os.path.exists(filename)
    df.to_csv(filename, mode='a', index=False, header=not file_exists)