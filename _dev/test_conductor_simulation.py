import conductor_simulation

if __name__ == "__main__":
    M, iterations = conductor_simulation.simulate_example_ECSE543_A1_Q3_monitored(w=1)
    M2, iterations2 = conductor_simulation.simulate_example_ECSE543_A1_Q3_monitored(use_gauss_seidel=True)
    print("Debug Point")