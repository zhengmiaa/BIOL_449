import numpy as np
import os
import pandas as pd
from PIL import Image

def calculate_overlap_area(ab_masks, cell_masks, pixel_to_micrometer):
    # Find the overlap region
    overlap_region = (ab_masks != 0) & (cell_masks != 0)
    
    # Count the number of pixels in the overlap
    overlap_pixels = np.count_nonzero(overlap_region)
    
    # Convert to physical area
    overlap_area_um2 = overlap_pixels * (pixel_to_micrometer ** 2)
    
    return overlap_area_um2

# parameters should be in um
def is_left(x, y, midpoint, perpendicular_vector):
        # vector of (x,y) relative to midpoint
        vec_x = x - midpoint[0]
        vec_y = y - midpoint[1]
        # cal cross product, decide direction
        cross = vec_x * perpendicular_vector[1] - vec_y * perpendicular_vector[0]
        # if cross > 0, point vector is on the left of perpendicular vector (counterclockwise from perp to point)
        # Q? using cross < 0 produces the correct Left Or Right? Why?
        return cross < 0

# Function to analyze the mask
def analyze_mask(cell_npy_path, ab_npy_path, pixel_to_micrometer, age, sub_midpoint, perpendicular_vector):
 # Load cell and ab masks
    cell_data = np.load(cell_npy_path, allow_pickle=True).item()
    ab_data = np.load(ab_npy_path, allow_pickle=True).item()

    # Extract masks
    cell_masks = np.array(cell_data['masks'])
    ab_masks = np.array(ab_data['masks'])

    intra_overlap_area_um2 = calculate_overlap_area(ab_masks, cell_masks, pixel_to_micrometer)
    print(f"Total intracellular ab accumulation area: {intra_overlap_area_um2:.2f} μm²")

    if (age == '6 weeks'):
        # remove any intracellular ab (overlap of ab and cell) for 6 week
        ab_masks[cell_masks != 0] = 0

    # -------- NEW!! stat about each side of SUB --------
    left_cell, right_cell = {'count':0, 'area':0.0, 'area_list':[], 'centroids':[]}, {'count':0, 'area':0.0, 'area_list':[], 'centroids':[]}
    left_ab, right_ab = {'count':0, 'area':0.0, 'area_list':[], 'centroids':[]}, {'count':0, 'area':0.0, 'area_list':[], 'centroids':[]}

    cell_mask_area_dist_um2 = []
    ab_mask_area_dist_um2 = []

    # STAT about CELL
    for cell_id in np.unique(cell_masks):
        if cell_id == 0:
            continue
            
        mask = cell_masks == cell_id
        area_pixels = np.sum(mask)
        area_um = area_pixels * (pixel_to_micrometer**2)

        # global cell area dist
        cell_mask_area_dist_um2.append(area_um)
        
        # cell centroid in pixel
        y_coords, x_coords = np.where(mask)
        centroid = (np.mean(x_coords), np.mean(y_coords))
        centroid_um = (centroid[0] * pixel_to_micrometer, centroid[1] * pixel_to_micrometer)

        # determine side
        if is_left(*centroid_um, sub_midpoint, perpendicular_vector):
            left_cell['count'] += 1
            left_cell['area'] += area_um
            left_cell['area_list'].append(area_um)
            left_cell['centroids'].append(centroid_um)
        else:
            right_cell['count'] += 1
            right_cell['area'] += area_um
            right_cell['area_list'].append(area_um)
            right_cell['centroids'].append(centroid_um)

    # STAT about AB
    for ab_id in np.unique(ab_masks):
        if ab_id == 0:
            continue
            
        mask = ab_masks == ab_id
        area_pixels = np.sum(mask)
        area_um = area_pixels * (pixel_to_micrometer**2)
        # global ab area distribution
        ab_mask_area_dist_um2.append(area_um)
        
        # AB centroid
        y_coords, x_coords = np.where(mask)
        centroid = (np.mean(x_coords), np.mean(y_coords))
        centroid_um = (centroid[0] * pixel_to_micrometer, centroid[1] * pixel_to_micrometer)
        
        # side
        if is_left(*centroid_um, sub_midpoint, perpendicular_vector):
            left_ab['count'] += 1
            left_ab['area'] += area_um
            left_ab['area_list'].append(area_um)
            left_ab['centroids'].append(centroid_um)
        else:
            right_ab['count'] += 1
            right_ab['area'] += area_um
            right_ab['area_list'].append(area_um)
            right_ab['centroids'].append(centroid_um)

    # Cal pixel overlapping
    overlap_mask = (ab_masks != 0) & (cell_masks != 0)
    overlap_pixels = np.sum(overlap_mask)
    y_coords, x_coords = np.where(overlap_mask)
    
    # in Pixels!!
    left_overlap = np.sum([is_left(x, y, sub_midpoint, perpendicular_vector) for x, y in zip(x_coords, y_coords)])
    right_overlap = overlap_pixels - left_overlap
    
    # 转换为实际面积
    left_overlap *= (pixel_to_micrometer**2)
    right_overlap *= (pixel_to_micrometer**2)

    print("\n=== Overall Stat ===")
    print(f"Cell Count: {left_cell['count'] + right_cell['count']}")
    print(f"AB Area: {left_ab['area'] + right_ab['area']:.2f} μm²")
    print(f"Intra AB Area: {intra_overlap_area_um2:.2f} μm²")

    print("\n=== Left Stat ===")
    print(f"Left Cell Count: {left_cell['count']}")
    print(f"Left Cell Area: {left_cell['area']:.2f} μm²")
    print(f"Left AB Count: {left_ab['count']}")
    print(f"Left AB Area: {left_ab['area']:.2f} μm²")

    print("\n=== Right Stat ===")
    print(f"Right Cell Count: {right_cell['count']}")
    print(f"Right Cell Area: {right_cell['area']:.2f} μm²")
    print(f"Right AB Count: {right_ab['count']}")
    print(f"Right AB Area: {right_ab['area']:.2f} μm²")

    return {
        'cell_global_dist': cell_mask_area_dist_um2,
        'ab_global_dist': ab_mask_area_dist_um2,
        'left': {
            'cell_dist': left_cell['area_list'],
            'ab_dist': left_ab['area_list'],
            'cell_centroids': left_cell['centroids'],
            'ab_centroids': left_ab['centroids'],
            'cell_count': left_cell['count'],
            'cell_area': left_cell['area'],
            'plaques': left_ab['count'],
            'plaque_area': left_ab['area']
        },
        'right': {
            'cell_dist': right_cell['area_list'],
            'ab_dist': right_ab['area_list'],
            'cell_centroids': right_cell['centroids'],
            'ab_centroids': right_ab['centroids'],
            'cell_count': right_cell['count'],
            'cell_area': right_cell['area'],
            'plaques': right_ab['count'],
            'plaque_area': right_ab['area']
        },
        'total': {
            'cell_count': left_cell['count'] + right_cell['count'],
            'cell_area': left_cell['area'] + right_cell['area'],
            'plaques': left_ab['count'] + right_ab['count'],
            'plaque_area': left_ab['area'] + right_ab['area'],
            'intra_plaque': intra_overlap_area_um2
        }
    }

    # # Calculate total area in pixels
    # cell_mask_area_pixels = np.count_nonzero(cell_masks)
    # ab_mask_area_pixels = np.count_nonzero(ab_masks)

    # # Normalize to micrometers
    # cell_mask_area_um2 = cell_mask_area_pixels * (pixel_to_micrometer ** 2)
    # print(f"Total cell mask area: {cell_mask_area_um2:.2f} μm²")
    # ab_mask_area_um2 = ab_mask_area_pixels * (pixel_to_micrometer ** 2)
    # print(f"Total ab mask area: {ab_mask_area_um2:.2f} μm²")

    # # Get unique mask values excluding zero
    # unique_cell_masks = np.unique(cell_masks)
    # unique_ab_masks = np.unique(ab_masks)

    # cell_mask_area_dist_um2 = []
    # ab_mask_area_dist_um2 = []

    # # Calculate individual mask areas in micrometers²
    # for mask_value in unique_cell_masks:
    #     if mask_value != 0:
    #         area_pixels = np.count_nonzero(cell_masks == mask_value)
    #         area_um2 = area_pixels * (pixel_to_micrometer ** 2)
    #         cell_mask_area_dist_um2.append(area_um2)

    # for mask_value in unique_ab_masks:
    #     if mask_value != 0:
    #         area_pixels = np.count_nonzero(ab_masks == mask_value)
    #         area_um2 = area_pixels * (pixel_to_micrometer ** 2)
    #         ab_mask_area_dist_um2.append(area_um2)

    # # Print results

    

    # return (
    #     cell_mask_area_um2, 
    #     len(cell_mask_area_dist_um2), 
    #     cell_mask_area_dist_um2, 
    #     ab_mask_area_um2, 
    #     len(ab_mask_area_dist_um2), 
    #     ab_mask_area_dist_um2,
    #     intra_overlap_area_um2
    # )



# # Main function
# if __name__ == "__main__":
#     # Get image paths from user input
#     cell_npy_path = input("Please enter the path to your CELL BODY segmentation mask image file (e.g., C:/path/to/image_seg.npy): ")
#     ab_npy_path = input("Please enter the path to your AMYLOID BETA segmentation mask image file (e.g., C:/path/to/image_seg.npy): ")
#     image_path_sub = input("Please enter the path to your subiculum mask image file (e.g., C:/path/to/image.tif): ")

#     # Check if files exist
#     if not (os.path.isfile(cell_npy_path) and os.path.isfile(ab_npy_path) and os.path.isfile(image_path_sub)):
#         print("Error: One or more specified files do not exist. Please check the paths and try again.")
#         exit(1)

#     # Get physical and pixel dimensions
#     try:
#         physical_width = float(input("Enter the physical width of the image in micrometers (e.g., 1286.15): "))
#         # physical_height = float(input("Enter the physical height of the image in micrometers (e.g., 810.14): "))
#         pixel_width = int(input("Enter the pixel width of the image (e.g., 7487): "))
#         # pixel_height = int(input("Enter the pixel height of the image (e.g., 4716): "))
#     except ValueError:
#         print("Error: Please enter valid numerical values for dimensions.")
#         exit(1)

#     # Calculate pixel-to-micrometer conversion factor
#     pixel_to_micrometer = physical_width / pixel_width  # Assumes square pixels
#     print(f"Calculated pixel size: {pixel_to_micrometer:.5f} μm/pixel")

#     # Calculate the total area of the masks
#     cell_mask_area, unique_cell_masks, cell_mask_area_dist, ab_mask_area, unique_ab_masks, ab_mask_area_dist = analyze_mask(
#         cell_npy_path, ab_npy_path, pixel_to_micrometer
#     )
#     sub_area = analyze_mask_sub(image_path_sub, pixel_to_micrometer)   # Assuming `analyze_mask_sub` handles its own scaling

#     # Save results
#     save_results(
#         cell_mask_area, unique_cell_masks, cell_mask_area_dist,
#         ab_mask_area, unique_ab_masks, ab_mask_area_dist, sub_area
#     )
