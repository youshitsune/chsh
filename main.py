import streamlit as st
import qiskit
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister, transpile, Aer
from qiskit.tools.monitor import job_monitor, backend_monitor, backend_overview
import matplotlib.pyplot as plt
import numpy as np

def make_chsh_circuit(theta_vec):
    chsh_circuits = []

    for theta in theta_vec:
        obs_vec = ['00', '01', '10', '11']
        for el in obs_vec:
            qc = QuantumCircuit(2, 2)
            qc.h(0)
            qc.cx(0,1)
            qc.ry(theta, 0)
            for a in range(2):
                if el[a] == '1':
                    qc.h(a)
            qc.measure(range(2), range(2))
            chsh_circuits.append(qc)
    return chsh_circuits

def compute_chsh_witness(counts):
    CHSH1 = []
    CHSH2 = []
    for i in range(0, len(counts), 4):
        theta_dict = counts[i:i + 4]
        zz = theta_dict[0]
        zx = theta_dict[1]
        xz = theta_dict[2]
        xx = theta_dict[3]

        no_shots = sum(xx[y] for y in xx)
        chsh1 = 0
        chsh2 = 0

        for element in zz:
            parity = (-1)**(int(element[0])+int(element[1]))
            chsh1+=parity*zz[element]
            chsh2+=parity*zz[element]
    
        for element in zx:
            parity = (-1)**(int(element[0])+int(element[1]))
            chsh1+= parity*zx[element]
            chsh2-= parity*zx[element]
    
        for element in xz:
            parity = (-1)**(int(element[0])+int(element[1]))
            chsh1-= parity*xz[element]
            chsh2+= parity*xz[element]
    
        for element in xx:
            parity = (-1)**(int(element[0])+int(element[1]))
            chsh1+= parity*xx[element]
            chsh2+= parity*xx[element]
    
        CHSH1.append(chsh1/no_shots)
        CHSH2.append(chsh2/no_shots)

    return CHSH1, CHSH2
def run():
    sim = Aer.get_backend("aer_simulator")
    number_of_thetas = 15
    theta_vec = np.linspace(0,2*np.pi,number_of_thetas)
    my_chsh_circuits = make_chsh_circuit(theta_vec)

    result_ideal = sim.run(my_chsh_circuits).result()

    CHSH1, CHSH2 = compute_chsh_witness(result_ideal.get_counts())
    
    plt.figure(figsize=(12,8))
    plt.plot(theta_vec, CHSH1,'o-',label = 'CHSH1 Noiseless')
    plt.plot(theta_vec, CHSH2, 'o-', label='CHSH2 Noiseless')

    st.set_option('deprecation.showPyplotGlobalUse', False)
    st.pyplot()

if st.button("Run simulation"):
    run()

