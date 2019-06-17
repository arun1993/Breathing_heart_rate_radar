from imp import load_source
WalabotAPI = load_source('WalabotAPI','/usr/share/walabot/python/WalabotAPI.py')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import datetime as dt
import numpy as np
# Parameters
x_len = 200         # Number of points to display
y_range = [0.0000001,0.000001]  # Range of possible Y values to display

# Create figure for plotting
# fig = plt.figure()
# ax = fig.add_subplot(1, 1, 1)
# xs = list(range(0, 200))
# ys = [0] * x_len
# ax.set_ylim(y_range)
#
# line, = ax.plot(xs, ys)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
xs = []
ys = []

# Add labels
plt.title('Breathing rate')
plt.xlabel('Samples')
plt.ylabel('Breathing energy')

wlbt = WalabotAPI
#wlbt.Clean()
wlbt.Init()
wlbt.SetSettingsFolder()
wlbt.ConnectAny()


wlbt.SetProfile(wlbt.PROF_SENSOR_NARROW)
wlbt.SetDynamicImageFilter(wlbt.FILTER_TYPE_DERIVATIVE)
wlbt.SetArenaTheta(-0.1, 0.1, 10)
wlbt.SetArenaPhi(-0.1, 0.1, 10)
wlbt.SetArenaR(20, 80, 0.2)
wlbt.Start()
wlbt.StartCalibration()
#wlbt.SetThreshold(100)

def animate(i, ys):

    wlbt.Trigger() # trigger walabot
    ener = wlbt.GetImageEnergy() # use image energy of walabot to detect breathing
    print(ener)
    # Add y to list
    ys.append(ener)

    # Limit y list to set number of items
    ys = ys[-x_len:]

    # Update line with new Y values
    line.set_ydata(ys)

    return line,

def animate2(i, xs, ys):
    xs.append(dt.datetime.now().strftime('%H:%M:%S.%f'))
    wlbt.Trigger()  # trigger walabot
    ener = wlbt.GetImageEnergy()  # use image energy of walabot to detect breathing
    print(ener)
    ys.append(ener)

    # Limit x and y lists to 20 items
    xs = xs[-20:]
    ys = ys[-20:]

    # Draw x and y lists
    ax.clear()
    ax.plot(xs, ys)

    # Format plot
    plt.xticks(rotation=45, ha='right')
    plt.subplots_adjust(bottom=0.30)
    plt.title('Time')
    plt.ylabel('Breathing energy')


# Set up plot to call animate() function periodically
ani = animation.FuncAnimation(fig, animate2, fargs=(xs, ys), interval=300)
plt.show()

# # Set up plot to call animate() function periodically
# ani = animation.FuncAnimation(fig,
#     animate,
#     fargs=(ys,),
#     interval=300,
#     blit=True)
# plt.show()

# pylive = load_source('pylive','C:/Users/Arun/PycharmProjects/walabot/pylive.py')
# from pylive import live_plotter
#
# import numpy as np
#
# size = 100
# x_vec = np.linspace(0,1,size+1)[0:-1]
# y_vec = np.random.randn(len(x_vec))
# line1 = []
# while True:
#     #wlbt.Trigger()# trigger walabot
#     #ener = wlbt.GetImageEnergy()# use image energy of walabot to detect breathing
#     #print(ener)
#     #rand_val = ener
#     y_vec[-1] = np.random.randn(1)
#     line1 = live_plotter(x_vec,y_vec,line1)
#     y_vec = np.append(y_vec[1:],0.0)
