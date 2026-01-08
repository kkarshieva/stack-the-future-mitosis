from image import smiles_to_svg

smiles = "Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1"

try:
    svg_bytes = smiles_to_svg(smiles)
    svg_str = svg_bytes.decode()  # convert bytes to string
    print("SVG string generated successfully!\n")
    print(svg_str[:500])  # print first 500 chars to check
except Exception as e:
    print("Error generating SVG:", e)

with open("test_molecule.svg", "wb") as f:
    f.write(svg_bytes)
print("SVG file 'test_molecule.svg' written successfully.")
