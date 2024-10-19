import finite_difference_mesh

if __name__ == "__main__":
    for i in range(2, 16):
        c, pos_node = finite_difference_mesh.generate_fd_analysis_circuit(i, to_csv=True)
        v = c.voltages()