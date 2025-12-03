#!/usr/bin/env python3
"""
ASCII Ramp Generator - Calculates pixel density for ASCII characters
and generates optimal brightness ramps for terminal rendering.
"""

from PIL import Image, ImageDraw, ImageFont


def calculate_char_density(char: str, font_size: int = 20) -> tuple[str, int, float]:
    """
    Render a character and count filled pixels.

    Returns tuple of (character, pixel_count, density_ratio)
    """
    # Create image with white background
    img_size = (font_size * 2, font_size * 2)
    img = Image.new('L', img_size, color=255)
    draw = ImageDraw.Draw(img)

    # Use default monospace font
    font = ImageFont.load_default()

    # Draw character in black
    draw.text(
        (font_size // 2, font_size // 2),
        char,
        fill=0,
        font=font,
    )

    # Count black pixels (filled)
    pixels = list(img.getdata())
    black_pixels = sum(1 for p in pixels if p < 128)
    total_pixels = len(pixels)
    density = black_pixels / total_pixels if total_pixels > 0 else 0.0

    return char, black_pixels, density


def generate_ramp(chars: str) -> str:
    """
    Generate ASCII ramp sorted by pixel density (dark to bright).
    """
    results = []

    for char in chars:
        if char.isprintable() and not char.isspace() or char == ' ':
            char_info, pixel_count, density = calculate_char_density(char)
            results.append((char_info, pixel_count, density))

    # Sort by density (ascending = dark to bright)
    results.sort(key=lambda x: x[1])

    # Print detailed results
    print("Character Pixel Density Analysis:")
    print("-" * 50)
    print(f"{'Char':<6} {'Pixels':<10} {'Density':<10}")
    print("-" * 50)

    for char, pixels, density in results:
        # Escape special characters for display
        display_char = repr(char) if char in [' ', '\t'] else char
        print(f"{display_char:<6} {pixels:<10} {density:<10.4f}")

    # Generate ramp string
    ramp = ''.join(char for char, _, _ in results)
    return ramp


def main():
    # Standard ASCII printable characters commonly used in ASCII art
    # Excluding letters/numbers to focus on symbols and punctuation
    standard_chars = " .`',-_:;^~\"<>!=+*?/\\|()[]{}#$%&@"

    print("Analyzing standard ASCII characters...")
    print()

    ramp = generate_ramp(standard_chars)

    print()
    print("=" * 50)
    print("Generated ASCII Ramp (dark to bright):")
    print("=" * 50)
    print(f'static char RAMP[] = "{ramp}";')
    print()
    print(f"Ramp length: {len(ramp)} characters")
    print()

    # Also generate a more conservative ramp with fewer characters
    conservative_chars = " .',-:;~=+*#%@"
    print("=" * 50)
    print("Conservative Ramp (classic style, more levels):")
    print("=" * 50)
    ramp_conservative = generate_ramp(conservative_chars)
    print()
    print(f'static char RAMP[] = "{ramp_conservative}";')
    print()


if __name__ == "__main__":
    main()
