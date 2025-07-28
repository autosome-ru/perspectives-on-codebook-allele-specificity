import os
import numpy as np

INPUT_DIR = 'motifs'
PPM_DIR = 'ppms'
PWM_DIR = 'pwms'

os.makedirs(PPM_DIR, exist_ok=True)
os.makedirs(PWM_DIR, exist_ok=True)

def parse_matrix_file(filepath):
    tf_name = None
    with open(filepath) as f:
        lines = f.readlines()
    for line in lines:
        if line.startswith('TF'):
            parts = line.strip().split()
            if len(parts) > 1:
                tf_name = parts[1]
    data_start = None
    for i, line in enumerate(lines):
        if line.startswith('Pos'):
            data_start = i + 1
            break
    if data_start is None:
        raise ValueError(f"File {filepath} is missing header line starting with 'Pos'")
    matrix_lines = [l for l in lines[data_start:] if l.strip() and l[0].isdigit()]
    matrix = [list(map(float, l.strip().split()[1:])) for l in matrix_lines]
    return tf_name, np.array(matrix)

def write_ppm(tf_name, matrix, output_path):
    with open(output_path, 'w') as out:
        out.write(f'>{tf_name}\n')
        for row in matrix:
            out.write('\t'.join(f'{x:.6f}' for x in row) + '\n')

def write_pwm(tf_name, matrix, output_path):
    with open(output_path, 'w') as out:
        out.write(f'>{tf_name}\n')
        pwm_matrix = np.log((matrix + 1e-5) / 0.25)
        for row in pwm_matrix:
            out.write('\t'.join(f'{x:.12f}' for x in row) + '\n')

def main():
    for fname in os.listdir(INPUT_DIR):
        if not fname.endswith('.txt'):
            continue
        in_path = os.path.join(INPUT_DIR, fname)
        tf_name_simple = fname.split('_')[0]
        base_out = os.path.splitext(tf_name_simple)[0]
        try:
            inner_name, matrix = parse_matrix_file(in_path)
        except Exception as e:
            print(f"Skipping {fname} due to error: {e}")
            continue
        tf_out_name = base_out
        ppm_out = os.path.join(PPM_DIR, f"{tf_out_name}.ppm")
        pwm_out = os.path.join(PWM_DIR, f"{tf_out_name}.pwm")
        write_ppm(tf_out_name, matrix, ppm_out)
        write_pwm(tf_out_name, matrix, pwm_out)

if __name__ == '__main__':
    main()

