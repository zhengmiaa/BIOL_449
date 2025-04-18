from PIL import Image
import numpy as np
from typing import Dict, Tuple, Optional

def calculate_end_connecting_line(
    L_end_pixel: Tuple[float, float], 
    R_end_pixel: Tuple[float, float],
    pixel_to_micrometer
) -> Dict[str, float]:
    # coordinates of the 2 ends of sub in um
    x1 = L_end_pixel[0] * pixel_to_micrometer
    y1 = L_end_pixel[1] * pixel_to_micrometer
    x2 = R_end_pixel[0] * pixel_to_micrometer
    y2 = R_end_pixel[1] * pixel_to_micrometer
    
    result = {
        'slope': None,
        'intercept': None, #um
        'direction_vector': (0.0, 0.0), #um
        'is_vertical': False,
        'is_horizontal': False,
        'midpoint_um': ((x1 + x2)/2, (y1 + y2)/2)
    }
    
    dx = x2 - x1
    dy = y2 - y1    
    # In case the line is vertical
    if abs(dx) < 1e-6: # do not use if dx == 0, bc dx is float maybe close to 0 but !=0
        return {
            **result,
            'is_vertical': True,
            'intercept': x1, #um
            'direction_vector': (0.0, 1.0)
        }
    
    slope = dy / dx
    intercept = y1 - slope * x1 #um
    
    return {
        **result,
        'slope': slope,
        'intercept': intercept,
        'direction_vector': (dx, dy),   #um
        'is_horizontal': abs(slope) < 1e-6
    }

def calculate_perpendicular_line(
    connect_line: Dict[str, float]
) -> Dict[str, float]:
    x_mid, y_mid = connect_line['midpoint_um']
    
    base_result = {
        'slope': None,
        'intercept': None,  #um
        'direction_vector': (0.0, 0.0), #um
        'is_vertical': False,
        'is_horizontal': False
    }
    
    # case that original line is vertical (should not be!)
    if connect_line['is_vertical']:
        return {
            **base_result,
            'slope': 0.0,
            'intercept': y_mid,
            'is_horizontal': True,
            'direction_vector': (1.0, 0.0)
        }
    
    # case that the original line is horizontal (unlikely)
    if connect_line['is_horizontal']:
        return {
            **base_result,
            'is_vertical': True,
            'intercept': x_mid,
            'direction_vector': (0.0, 1.0)
        }
    
    # General cases (neither vertical nor horizontal)
    if connect_line['slope'] is not None:
        perpendicular_slope = -1 / connect_line['slope']
    else:
        perpendicular_slope = None  # Should not happen, but added for safety

    perpendicular_intercept = y_mid - perpendicular_slope * x_mid
    dx, dy = connect_line['direction_vector']
    perpendicular_vector = (-dy, dx)  # rotate 90 degree
    
    return {
        **base_result,
        'slope': perpendicular_slope,
        'intercept': perpendicular_intercept,
        'direction_vector': perpendicular_vector
    }

def analyze_mask_sub(
    image_path_sub: str,
    pixel_to_micrometer: float,
    left_sub_end: Tuple[float, float],
    right_sub_end: Tuple[float, float]
) -> Dict:
    """
    return data structure
    {
        "geometry": {
            "area_um2": ...,
            "width_um": ...,
            "height_um": ...,
            "bounding_box": {
                "min_row": ..., "max_row": ...,
                "min_col": ..., "max_col": ...
            }
        },
        "midline": {
            "slope": ...,
            "intercept": ...,
            "direction_vector": (...),
            "is_vertical": ...,
            "is_horizontal": ...,
            'midpoint_um': ...
        },
        "dividing_line": {
            "slope": ...,
            "intercept": ...,
            "direction_vector": (...),
            "is_vertical": ...,
            "is_horizontal": ...
        }
    }
    """
    # load image
    with Image.open(image_path_sub) as img:
        img_array = np.array(img)
    
    binary_mask = (img_array == 255).astype(np.uint8)
    white_coords = np.argwhere(binary_mask == 1)
    
    if white_coords.size == 0:
        raise ValueError("No white region found in the image")
    
    # calculate the geometry
    min_row, min_col = white_coords.min(axis=0)
    max_row, max_col = white_coords.max(axis=0)
    
    geometry = {
        "area_um2": binary_mask.sum() * (pixel_to_micrometer ** 2),
        "width_um": (max_col - min_col + 1) * pixel_to_micrometer,
        "height_um": (max_row - min_row + 1) * pixel_to_micrometer,
        "bounding_box": {
            "min_row": min_row,
            "max_row": max_row,
            "min_col": min_col,
            "max_col": max_col
        }
    }
    
    # cal midline (connecting two ends) in um!!
    connect_line = calculate_end_connecting_line(left_sub_end, right_sub_end, pixel_to_micrometer)
    
    # cal dividing line
    perpendicular_line = calculate_perpendicular_line(connect_line)

    return {
        "geometry": geometry,
        "midline": connect_line,
        "dividing_line": perpendicular_line
    }

