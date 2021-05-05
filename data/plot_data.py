import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd


data = pd.read_csv("data_motor400.csv", sep=';')
nSamples = 100
#iteration [];theta [°];theta [rad];Couple [Nm];omega [rad/s]
iteration = data["iteration []"].values[1:nSamples]
theta_degree = data["theta [°]"].values[1:nSamples]
theta_rad = data["theta [rad]"].values[1:nSamples]
torque = data["Couple [Nm]"].values[1:nSamples]*1000
omega = data["omega [rad/s]"].values[1:nSamples]




#filename = "Couple-Omega (motor400).txt"
#X = np.loadtxt(filename)
# t = np.arange(len(X))
# fig, ax = plt.subplots(1,2)
# ax[0].plot(omega,torque)
# ax[0].set_title("Couple en fonction de omega")
# ax[1].plot(iteration,torque)
# ax[1].set_title("Couple à chaque itération")
# plt.figure();
#plt.plot(omega*180/np.pi,torque);


# Angle [°] à chaque itération
fig_angle_degree_iteration = plt.figure()
ax = fig_angle_degree_iteration.add_subplot(1, 1, 1)
ax.plot(iteration, theta_degree, color='tab:blue')
ax.set_title('Angle de rotation du rotor à chaque itération')
ax.legend(['Position angulaire [°]'])
ax.set_xlabel('Iteration []')
ax.set_ylabel('Position angulaire' [°]')
plt.show()

# Angle [rad] à chaque itération
fig_angle_rad_iteration = plt.figure()
ax = fig_angle_rad_iteration.add_subplot(1, 1, 1)
ax.plot(iteration, theta_rad, color='tab:blue')
ax.set_title('Angle de rotation du rotor à chaque itération')
ax.legend(['Position angulaire [rad]'])
ax.set_xlabel('Iteration []')
ax.set_ylabel('Position angulaire [rad]')
plt.show()

# Vitesse angulaire [rad/s] à chaque itération
fig_omega_iteration = plt.figure()
ax = fig_omega_iteration.add_subplot(1, 1, 1)
ax.plot(iteration, omega, color='tab:blue')
ax.set_title('Vitesse angulaire du rotor à chaque itération')
ax.legend(['Vitesse angulaire [rad/s]'])
ax.set_xlabel('Iteration []')
ax.set_ylabel('Omega [rad/s]')
plt.show()

# Couple à chaque itération
avTorque = np.mean(torque)
averageTorque = np.zeros(nSamples-1);
for i in range(nSamples-1):
    averageTorque[i] = avTorque
fig_torque_iteration = plt.figure()
ax = fig_torque_iteration.add_subplot(1, 1, 1)
ax.plot(iteration, torque, color='tab:blue')
ax.plot(iteration, averageTorque, color='tab:orange')
ax.set_title('Couple à chaque itération')
ax.legend(['Couple', 'Couple moyen'])
ax.set_xlabel('Iteration []')
ax.set_ylabel('Couple [mNm]')
plt.show()

#Puissance
averagePower = np.zeros(nSamples-1);
sumPower = 0;
for i in range(nSamples-1):
    sumPower += torque[i] * omega[i]
    averagePower[i] = sumPower/i
    
fig_power = plt.figure()
ax = fig_power.add_subplot(1, 1, 1)
ax.plot(iteration, torque*omega, color='tab:cyan')
ax.plot(iteration, averagePower, color='tab:orange')
ax.set_title('Puissance à chaque itération')
ax.legend(['Puissance', 'Puissance moyenne'])
ax.set_xlabel('Iteration []')
ax.set_ylabel('Puissance [mW]')
plt.show()

# Couple en fonction de la vitesse angulaire
fig_torque_angular_velocity = plt.figure()
ax = fig_torque_angular_velocity.add_subplot(1, 1, 1)
ax.plot(omega, torque, color='tab:blue')
ax.plot(omega, averageTorque, color='tab:orange')
ax.set_title('Couple en fonction de la vitesse angulaire')
ax.legend(['Couple', 'Couple moyen'])
ax.set_xlabel('Omega [rad/s]')
ax.set_ylabel('Couple [mNm]')
plt.show()