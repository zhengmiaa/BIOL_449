import os
from cal_area import analyze_mask
from analyze_SUB import analyze_mask_sub
from save_results import save_results

# Main function
if __name__ == "__main__":
    # Get image paths from user input
    cell_npy_path = input("Please enter the path to your CELL BODY segmentation mask image file (e.g., C:/path/to/image_seg.npy): ").strip('"')
    ab_npy_path = input("Please enter the path to your AMYLOID BETA segmentation mask image file (e.g., C:/path/to/image_seg.npy): ").strip('"')
    sub_path = input("Please enter the path to your subiculum mask image file (e.g., C:/path/to/image.tif): ").strip('"')
    
    # Check if files exist
    if not (os.path.isfile(cell_npy_path) and os.path.isfile(ab_npy_path) and os.path.isfile(sub_path)):
        print("Error: One or more specified files do not exist. Please check the paths and try again.")
        exit(1)

    # Get physical and pixel dimensions
    try:
        physical_width = float(input("Enter the physical width of the image in micrometers (e.g., 1286.15): "))
        physical_height = float(input("Enter the physical height of the image in micrometers (e.g., 810.14): "))
        pixel_width = int(input("Enter the pixel width of the image (e.g., 7487): "))
        # pixel_height = int(input("Enter the pixel height of the image (e.g., 4716): "))
    except ValueError:
        print("Error: Please enter valid numerical values for dimensions.")
        exit(1)

    # Calculate pixel-to-micrometer conversion factor
    pixel_to_micrometer = physical_width / pixel_width  # Assumes square pixels
    print(f"Calculated pixel size: {pixel_to_micrometer:.5f} μm/pixel")

    left_sub_end_pixel = tuple(map(float, input("please enter the left end of the subiculum (e.g. 1000, 2000) in pixel: ").split(',')))
    right_sub_end_pixel = tuple(map(float, input("please enter the right end of the subiculum (e.g. 1000, 2000) in pixel: ").split(',')))

    age_choice = input("6-week-old segmentation is handled differently by removing intracellular ab accumulation. Choose Age (1. 6 weeks, 2. 10 weeks, 3. 6 months):")
    age_mapping = {'1': '6 weeks','2': '10 weeks','3': '6 months'}
    age = age_mapping.get(age_choice, "NA")

    # sub_result structure:
    # {
    #     "geometry": {
    #         "area_um2": ...,
    #         "width_um": ...,
    #         "height_um": ...,
    #         "bounding_box": {
    #             "min_row": ..., "max_row": ...,
    #             "min_col": ..., "max_col": ...
    #         }
    #     },
    #     "midline": {
    #         "slope": ...,
    #         "intercept": ...,
    #         "direction_vector": (...),
    #         "is_vertical": ...,
    #         "is_horizontal": ...,
    #         'midpoint_um': ...
    #     },
    #     "dividing_line": {
    #         "slope": ...,
    #         "intercept": ...,
    #         "direction_vector": (...),
    #         "is_vertical": ...,
    #         "is_horizontal": ...
    #     }
    # }
    
    sub_result = analyze_mask_sub(sub_path, pixel_to_micrometer, left_sub_end_pixel, right_sub_end_pixel)

    # mask result is of the following structure:
    #     {
    #     'cell_global_dist': [],
    #     'ab_global_dist': [],
    #     'left': {
    #         'cell_dist': [],
    #         'ab_dist': [],
    #         ...else...
    #     },
    #     'right': { ... },
    #     'total': { ... }
    # }
    mask_result = analyze_mask(
        cell_npy_path, ab_npy_path, pixel_to_micrometer, age, 
        sub_result["midline"]["midpoint_um"], sub_result["dividing_line"]["direction_vector"]
    )

    mouse_id = input("Enter Mouse ID: ")
    sex_choice = input("Choose Sex (1. Female, 2.Male): ")
    sex_mapping = {'1': 'Female', '2':'Male'}
    sex = sex_mapping.get(sex_choice, "NA")

    genotype_choice = input("Choose Genotype (1.WT, 2.5XFAD): ")
    genotype_mapping = {'1': 'WT', '2':'5XFAD'}
    genotype = genotype_mapping.get(genotype_choice, "NA")

    save_results(
        mouse_id, sex, genotype, age,
        pixel_to_micrometer, physical_width, physical_height,
        mask_result, 
        sub_result        
    )