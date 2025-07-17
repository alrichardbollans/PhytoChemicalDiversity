from PIL import Image, ImageDraw, ImageFont


def combine_images_side_by_side(
        image1_path,
        image2_path,
        output_path="combined_image.jpg",
        separation=10,
        text1=None,
        text2=None,
        font_path=None,
        font_size=20,
        side_by_side=True
):
    # Open the two images
    img1 = Image.open(image1_path)
    img2 = Image.open(image2_path)

    # Get dimensions for each image
    width1, height1 = img1.size
    width2, height2 = img2.size

    # Prepare for text height if any text is provided
    text_height = 0
    if text1 or text2:
        if font_path:
            font = ImageFont.truetype(font_path, font_size)
        else:
            font = ImageFont.load_default(font_size)

        # Estimate the height of the text
        text_height = font.getbbox("A")[3] + 10  # Add padding below text


    # Set the height to the maximum height of the two images + text height if applicable
    if side_by_side:
        combined_height = max(height1, height2) + text_height
        combined_width = width1 + width2 + separation

    else:
        combined_height = height1 + height2 + 2*text_height + separation
        combined_width = max(width1, width2)

    # Create a new blank image with the combined dimensions
    combined_image = Image.new("RGB", (combined_width, combined_height), color=(255, 255, 255))

    # Paste each image onto the new image
    combined_image.paste(img1, (0, 0))
    if side_by_side:
        combined_image.paste(img2, (width1+separation, 0))
    else:
        combined_image.paste(img2, (0, height1+ separation))



    # Draw the text under each image if provided
    draw = ImageDraw.Draw(combined_image)
    if text1:
        text_width = font.getbbox(text1)[2]  # Get the width of the text
        text_x = width1 // 2 - text_width // 2
        draw.text((text_x, height1 + 5), text1, fill="black", font=font)
    if text2:
        text_width = font.getbbox(text2)[2]  # Get the width of the text
        if side_by_side:
            text_x = width1 + separation + (width2 // 2) - text_width // 2
            draw.text((text_x, height2 + 5), text2, fill="black", font=font)
        else:
            text_x = width1 // 2 - text_width // 2
            draw.text((text_x, height1 + height2 + separation + 5), text2, fill="black", font=font)

    # Save the combined image
    combined_image.save(output_path)
    print(f"Combined image saved as {output_path}")

if __name__ == '__main__':

    # Example usage:
    combine_images_side_by_side(
        "outputs/species_in_study_native_dist.jpg",
        "outputs/species_in_families_native_dist.jpg",
        output_path="outputs/species_plots.jpg",
        separation=160,
        text1="(a) Species in Study",
        text2="(b) Species in Gentianales",
        # font_path="arial.ttf",  # Optional: specify your font path or leave as None for default
        font_size=160,
        side_by_side=False
    )
