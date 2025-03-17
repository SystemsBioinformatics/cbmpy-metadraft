import os
import math
from PIL import Image, ImageDraw, ImageFont

# Define a 256 RGB colour table (example with a gradient approach for simplicity)
def generate_rgb_colour_table():
    colour_table = []
    for i in range(256):
        # We are generating a gradient from black to red as an example
        r = i
        g = (255 - i) // 2
        b = (255 - i)
        colour_table.append((r, g, b))
    return colour_table


def extended_rgb_colour_table():
    RGB_COLOURS_136_= (
        (128,0,0),
        #(139,0,0),
        (165,42,42),
        (178,34,34),
        (220,20,60),
        (255,0,0),
        (255,99,71),
        (255,127,80),
        (205,92,92),
        (240,128,128),
        (233,150,122),
        (250,128,114),
        (255,160,122),
        (255,69,0),
        (255,140,0),
        (255,165,0),
        (255,215,0),
        (184,134,11),
        (218,165,32),
        (238,232,170),
        (189,183,107),
        (240,230,140),
        (128,128,0),
        (255,255,0),
        (154,205,50),
        (85,107,47),
        (107,142,35),
        (124,252,0),
        (127,255,0),
        (173,255,47),
        (0,100,0),
        (0,128,0),
        (34,139,34),
        (0,255,0),
        (50,205,50),
        (144,238,144),
        (152,251,152),
        (143,188,143),
        (0,250,154),
        (0,255,127),
        (46,139,87),
        (102,205,170),
        (60,179,113),
        (32,178,170),
        (47,79,79),
        (0,128,128),
        (0,139,139),
        (0,255,255),
        (0,255,255),
        (224,255,255),
        (0,206,209),
        (64,224,208),
        (72,209,204),
        (175,238,238),
        (127,255,212),
        (176,224,230),
        (95,158,160),
        (70,130,180),
        (100,149,237),
        (0,191,255),
        (30,144,255),
        (173,216,230),
        (135,206,235),
        (135,206,250),
        (25,25,112),
        (0,0,128),
        (0,0,139),
        (0,0,205),
        (0,0,255),
        (65,105,225),
        (138,43,226),
        (75,0,130),
        (72,61,139),
        (106,90,205),
        (123,104,238),
        (147,112,219),
        (139,0,139),
        (148,0,211),
        (153,50,204),
        (186,85,211),
        (128,0,128),
        (216,191,216),
        (221,160,221),
        (238,130,238),
        (255,0,255),
        (218,112,214),
        (199,21,133),
        (219,112,147),
        (255,20,147),
        (255,105,180),
        (255,182,193),
        (255,192,203),
        (250,235,215),
        (245,245,220),
        (255,228,196),
        (255,235,205),
        (245,222,179),
        (255,248,220),
        (255,250,205),
        (250,250,210),
        (255,255,224),
        (139,69,19),
        (160,82,45),
        (210,105,30),
        (205,133,63),
        (244,164,96),
        (222,184,135),
        (210,180,140),
        (188,143,143),
        (255,228,181),
        (255,222,173),
        (255,218,185),
        (255,228,225),
        (255,240,245),
        (250,240,230),
        (253,245,230),
        (255,239,213),
        (255,245,238),
        (245,255,250),
        (112,128,144),
        (119,136,153),
        (176,196,222),
        (230,230,250),
        (255,250,240),
        (240,248,255),
        (248,248,255),
        (240,255,240),
        (255,255,240),
        (240,255,255),
        (255,250,250),
        (105,105,105),
        (128,128,128),
        (169,169,169),
        (192,192,192),
        (211,211,211),
        (220,220,220),
    )
    return RGB_COLOURS_136_
    

# Split colours into light and dark
def split_light_dark(colour_table, method='luminance_opt1'):
    light_colours, dark_colours = [], []
    for colour in colour_table:
        # Calculate the brightness of the colour
        # Inspired by https://stackoverflow.com/questions/596216/formula-to-determine-perceived-brightness-of-rgb-color
        
        if method == 'luminance':
            # print('Using method:', 'luminance')
            avg_brightness = math.floor(0.2126*colour[0] + 0.7152*colour[1] + 0.0722*colour[2]) 
            threshold = 160
        elif method == 'luminance_opt1':
            # print('Using method:', 'luminance_opt1')
            avg_brightness = math.floor(0.299*colour[0] + 0.587*colour[1] + 0.144 * colour[2])
            threshold = 160
        elif method == 'luminance_opt2':
            # print('Using method:', 'luminance_opt2')
            avg_brightness = math.floor(math.sqrt(0.299*colour[0]**2 + 0.587*colour[1]**2 + 0.144*colour[2]**2))
            threshold = 160
        else: # fallback raw RGB average
            # print('Using method:', 'RGB avg')
            avg_brightness = sum(colour) // 3
            threshold = 200

        if avg_brightness > threshold:
            light_colours.append(colour)
        else:
            dark_colours.append(colour)
    return light_colours, dark_colours

# Create an image for each colour list
def create_colour_table_image(colours, filename):
    cell_size = 50
    width, height = cell_size * 4, cell_size * len(colours) + 100
    image = Image.new('RGB', (width, height), (255, 255, 255))
    draw = ImageDraw.Draw(image)

    # Load a font
    try:
        font = ImageFont.truetype("arial.ttf", 15)
    except IOError:
        font = ImageFont.load_default()

    for idx, colour in enumerate(colours):
        # Draw each colour block
        y0 = idx * cell_size
        draw.rectangle([0, y0, cell_size, y0 + cell_size], fill=colour, outline=(0,0,0))

        # Add the RGB text overlay
        text = f"{colour}"
        text_position = (cell_size * 1.5, y0 + cell_size // 2 - 10)
        draw.text(text_position, text, fill=(0, 0, 0), font=font)

    image.save(filename)

# Main execution
def main():
    
    #colour_table = generate_rgb_colour_table()
    colour_table = extended_rgb_colour_table()
    
    luminance_method = 'luminance_opt1'
    print('Using perceived brightness method:', luminance_method)
    light_colours, dark_colours = split_light_dark(colour_table, method=luminance_method)

    # Ensuring output directory exists
    #output_dir = 'output_images'
    #os.makedirs(output_dir, exist_ok=True)
    output_dir = os.getcwd()

    # Save images for each colour category
    create_colour_table_image(light_colours, os.path.join(output_dir, 'light_colours.png'))
    create_colour_table_image(dark_colours, os.path.join(output_dir, 'dark_colours.png'))
    
    with open('rgb_tables.py', 'w') as F:
        F.write("RGB_EXTENDED_COLOUR_TABLE = (\n")
        for r in colour_table:
            F.write('    {},\n'.format(r))
        F.write(')\n\n')
        F.write("RGB_LIGHT_COLOUR_TABLE = (\n")
        for r in light_colours:
            F.write('    {},\n'.format(r))
        F.write(')\n\n')
        F.write("RGB_DARK_COLOUR_TABLE = (\n")
        for r in dark_colours:
            F.write('    {},\n'.format(r))
        F.write(')\n\n')

        F.close()


if __name__ == "__main__":
    main()

