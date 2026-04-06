import sys
import csv
import numpy as np
import matplotlib.pyplot as plt

def plot_fields(csv_file, nv):
    # Lê o CSV
    x_vals = []
    y_vals = []
    p_vals = []

    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        next(reader)  # Pula o cabeçalho
        for row in reader:
            x_vals.append(float(row[0]))
            y_vals.append(float(row[1]))
            p_vals.append(float(row[4]))  # p é a 5ª coluna (índice 4)

    num_data = len(x_vals)
    print(f"Dados lidos: {num_data} pontos")
    print(f"Esperado: {nv * nv} pontos (nv={nv})")
    
    if num_data != nv * nv:
        print(f"Aviso: Número de pontos ({num_data}) não corresponde a nv²={nv*nv}")
        # Tenta descobrir o tamanho real da malha
        import math
        real_nv = int(math.sqrt(num_data))
        if real_nv * real_nv == num_data:
            print(f"Corrigindo nv para {real_nv}")
            nv = real_nv
        else:
            raise ValueError(f"Número de pontos {num_data} não é um quadrado perfeito")

    # Reorganiza em arrays 2D (nv x nv)
    x = np.array(x_vals).reshape((nv, nv))
    y = np.array(y_vals).reshape((nv, nv))
    p = np.array(p_vals).reshape((nv, nv))

    # Calcula dx e dy (assumindo grade uniforme)
    if nv > 1:
        dx = x[0, 1] - x[0, 0]
        dy = y[1, 0] - y[0, 0]
    else:
        dx = 1.0  # Valor padrão se nv=1
        dy = 1.0

    # Plota o campo de pressão
    fig, ax = plt.subplots(figsize=(8, 6))
    c = ax.pcolormesh(x, y, p, shading='auto', cmap='viridis')
    plt.colorbar(c, ax=ax, label='Pressão')

    # Adiciona linhas da malha (pretas e finas)
    for i in range(nv + 1):
        ax.axvline(x=i * dx, color='black', linewidth=0.5)
    for j in range(nv + 1):
        ax.axhline(y=j * dy, color='black', linewidth=0.5)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Campo de Pressão')
    ax.set_aspect('equal')  # Mantém proporção quadrada

    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Uso: python plot_pressure.py <arquivo_csv> <nv>")
        sys.exit(1)
    plot_fields(sys.argv[1], int(sys.argv[2]))