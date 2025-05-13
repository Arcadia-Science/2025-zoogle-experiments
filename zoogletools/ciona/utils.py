def hex_to_plotly_rgba(hex_color, alpha=1.0):
    """Convert a hex color to an RGBA string.

    Args:
        hex_color (str): The hex color code (e.g., "#FF0000" for red)
        alpha (float): The opacity value between 0 and 1

    Returns:
        str: RGBA color string for Plotly
    """
    hex_color = hex_color.lstrip("#")
    rgb = tuple(int(hex_color[i : i + 2], 16) for i in (0, 2, 4))
    return f"rgba({rgb[0]},{rgb[1]},{rgb[2]},{alpha})"
