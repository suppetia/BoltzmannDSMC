import cv2
import numpy as np


def process_image(filename, export_path):
    image = cv2.imread(filename)

    # convert to grayscale
    gray_image = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)


    mat = np.array(gray_image)

    # print(mat)

    mat_mask = mat < 255
    # print(mat_mask)

    np.savetxt(export_path, mat_mask.astype("int"), fmt="%1d",
               header=f"{mat_mask.shape[0]} {mat_mask.shape[1]}", comments="")
    print(f"converted '{filename}' to correct format and exported it to '{export_path}'")


if __name__ == "__main__":
    image_path = 'data/circle_1024.png'
    process_image(image_path, "matrix3.txt")
    image_path = 'data/circle_256.png'
    process_image(image_path, "matrix4.txt")
    image_path = 'data/circle_64_128.png'
    process_image(image_path, "matrix5.txt")
    image_path = 'data/circle_32_64.png'
    process_image(image_path, "matrix6.txt")



