import streamlit as st
import qiskit
from qiskit import QuantumCircuit, QuantumRegister, transpile, Aer
from qiskit_ibm_provider import IBMProvider
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

def run_real():
    number_of_thetas = 15
    theta_vec = np.linspace(0,2*np.pi,number_of_thetas)
    my_chsh_circuits = make_chsh_circuit(theta_vec)

    provider = IBMProvider(st.secrets["key"])
    backend = provider.get_backend("ibmq_lima")
    transpiled_circuits = transpile(my_chsh_circuits, backend)
    job_real = backend.run(transpiled_circuits, shots=8192)
    result_real = job_real.results()

    CHSH1, CHSH2 = compute_chsh_witness(result_real.get_counts())

    plt.figure(figsize=(12, 8))
    plt.plot(theta_vec, CHSH1, 'o-', label = 'CHSH1 Lima')
    plt.plot(theta_vec, CHSH2, 'o-', label = 'CHSH2 Lima')

    st.pyplot()

def run_simulation():
    sim = Aer.get_backend("aer_simulator")
    number_of_thetas = 15
    theta_vec = np.linspace(0,2*np.pi,number_of_thetas)
    my_chsh_circuits = make_chsh_circuit(theta_vec)

    result_ideal = sim.run(my_chsh_circuits).result()

    CHSH1, CHSH2 = compute_chsh_witness(result_ideal.get_counts())
    
    plt.figure(figsize=(12,8))
    plt.plot(theta_vec, CHSH1,'o-',label = 'CHSH1 Noiseless')
    plt.plot(theta_vec, CHSH2, 'o-', label='CHSH2 Noiseless')

    st.pyplot()

st.set_option('deprecation.showPyplotGlobalUse', False)
st.write("## This is showcase of CHSH")
st.write("It basically showcases that realism in physics is wrong. If it was correct than max number possible is 2, but here we are getting approx. 2.8.")
st.write("Source code: https://github.com/youshitsune/chsh")
sim, real = st.columns(2)
with sim:
    if st.button("Run simulation"):
        run_simulation()

with real:
    if st.button("Run on real quantum computer"):
        run_real()

    st.write("Running on real quantum computer can take a while!")


