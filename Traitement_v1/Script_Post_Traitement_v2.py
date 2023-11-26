import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import spectrogram


def average_elements(input_list, num_elements):
    res = []
    longr = len(input_list)

    for i in range(0, longr, num_elements):
        output_list = []
        if i+num_elements>=longr:
            num_elements = longr - i
        chunk = input_list[i:i+num_elements]
        average = sum(chunk) / num_elements 
        output_list.append(average)
        output_list = output_list * num_elements
        res = res + output_list

    return res

def sliding_average(input_list, num_elements):
    res = []
    longr = len(input_list)

    for i in range(0, longr):
        output_list = []
        if i+num_elements>=longr:
            num_elements = longr - i
        chunk = input_list[i:i+num_elements]
        average = sum(chunk) / num_elements 
        output_list.append(average)
        #output_list = output_list * num_elements
        res = res + output_list

    return res

def sliding_sum(input_list, num_elements):
    res = [0]*(num_elements)
    longr = len(input_list)

    for i in range(num_elements, longr):
        output_list = []
        # if i+num_elements>=longr:
        #     num_elements = longr - i
        chunk = input_list[i-num_elements:i]
        average = sum(chunk)  
        output_list.append(average)
        #output_list = output_list * num_elements
        res = res + output_list

    return res

def sum_elements(input_list, num_elements):
    res = []
    longr = len(input_list)
    start_index = num_elements

    for i in range(0, longr, num_elements):
        output_list = []
        if i+num_elements>=longr:
            num_elements = longr - i
        chunk = input_list[i:i+num_elements]
        average = sum(chunk) 
        output_list.append(average)
        output_list = output_list * num_elements
        res = res + output_list

    return res

def calculate_derivative(x, y):
    n = len(x)
    derivative = [0]

    for i in range(1, n - 1):
        dx = x[i + 1] - x[i - 1]
        dy = y[i + 1] - y[i - 1]
        derivative_value = dy / dx
        derivative.append(derivative_value)

    return derivative

def local_derivative(list_x,list_y,num_elements):
    res = []
    longr = len(list_x)
    for i in range(0, longr, num_elements):
        output_list = []
        if i+num_elements>=longr:
            num_elements = longr - i
        chunk_x = list_x[i:i+num_elements]
        chunk_y = list_y[i:i+num_elements]
        output_list = calculate_derivative(chunk_x,chunk_y)
        res = res + output_list

    return res



def plot_data(csv_file, csv_file2, sample_Period,delta):
    #sample_Period in ms
    y1 = [] # Head A (Voie 1)
    y2 = [] # Head B (Voie 2)
    y2_nh = [] # Head B (Voie 2) with head distance differen

    # ======================= Data Extraction =================================
    with open(csv_file, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header row if present

        # Extraction of data from specified file
        for row in reader:
            y1.append(float(row[0]))  # Extraction of the data from Head A
            # y2.append(float(row[1]))  # Extraction of the data from Head B

    with open(csv_file_2, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header row if present

        # Extraction of data from specified file
        for row in reader:
            y2.append(float(row[0]))  # Extraction of the data from Head B
            y2_nh.append(float(row[0])-85)        


    # ============================= Data Treatment ============================
    x = (np.linspace(0, sample_Period, len(y1)))*len(y1)/1000

    y_avrg = []
    theta = []
    d_theta = [0]
    dy1 = [0]
    dy2 = [0]
    dy_avrg = [0]



    # first extraction
    for i in range(0,len(y1)):
        y_avrg.append( (y1[i]+y2[i])/2 ) # average position put in a list 
        theta.append(  ((y2[i]-85-y1[i])*(10**-3)) /delta )  # !!!!! THE COEFFECIENT MIGHT CHANGE DEPENDING ON HOW !!!!!
                                                                                         # !!!!!      ANGLE AND SPEED ARE CALCULATED           !!!!!

    # Débruitage de théta ?
    theta = average_elements(theta, 10) 
    N = len(theta)
    theta_f = np.fft.fft(theta)
    Theta_A = np.abs(theta_f) / N

    #d_t_theta = local_derivative(x,theta,10)

    #print("fft = ",fft_theta)
        
    min_diff = min(theta)
    max_diff = max(theta)

    angle = 0
    rot = [0]
    for i in range(1,len(y1)):
        dy1.append((y1[i]-y1[i-1])/(sample_Period*(10**-3)))
        dy2.append((y2[i]-y2[i-1])/(sample_Period*(10**-3)))
        dy_avrg.append((y_avrg[i]-y_avrg[i-1])/(sample_Period*(10**-3)))
        d_theta.append(abs(10*(theta[i]-theta[i-1])  *18))
        angle += abs(10*(theta[i]-theta[i-1])  *18)
        rot.append(angle)

    rot_s = sliding_sum(d_theta,1000)
    rot_s = average_elements(rot_s,10)
    # for i in theta:
    #     rot.append(i*(360/(max_diff - min_diff)))


    #cleaned
    dy_avrg2 = average_elements(dy_avrg,10)
    d_theta2 = average_elements(d_theta,10)
    dy1 = average_elements(dy1,10)
    dy2 = average_elements(dy2,10)

    Fs = 1/(sample_Period*(10**-3)) # Sample frequency

    # ========================== Plotting of the different graphs ===================
    fig_1, ax1 = plt.subplots()

    #print("len of H_B = ", len(y1))
    ax1.plot(x,y_avrg,'-b')
    ax1.set_xlabel('t (s)')
    ax1.set_ylabel('Position mm', color='b')
    #plt.ylabel('Déplacement en mm')

    ax2 = ax1.twinx()

    ax2.plot(x, dy_avrg2,'-r')
    ax2.set_ylabel('Vitesse translation (mm/s)', color='r')


    #plt.grid(True)
    plt.title('Position et vitesse de translation en fonction du temps')

    fig_2, ax3 = plt.subplots()
    ax3.plot(x,theta,'-b')
    ax3.set_xlabel('t (s)')
    ax3.set_ylabel('theta (°)', color='b')

    # ax4= ax3.twinx()
    # ax4.plot(x,rot_s,'-r')
    # ax4.set_ylabel('rotation (°/s)',color='r')
    plt.title('Rotation et vitesse de rotation en fonction du temps')

    fig_3, ax5 = plt.subplots()
    ax5.plot(x,y1,'-b')
    ax5.set_xlabel('t (s)')
    ax5.set_ylabel('Position mm (Head A)', color='b')

    ax6= ax5.twinx()
    ax6.plot(x,y2,'-r')
    ax6.set_ylabel('Position mm (HeadB)',color='r')
    plt.title('Position des deux lasers en fonction du temps')

    fig_4, ax7 = plt.subplots()
    ax7.plot(x,rot,'-b')
    ax7.set_xlabel('t (s)')
    ax7.set_ylabel('Rotation totale (°) ?', color='b')

    ax8= ax7.twinx()
    ax8.plot(x,rot_s,'-r')
    ax8.set_ylabel('Vitesse de rotation ? (°/s)',color='r')
    plt.title('° totals et °/s ')

    print("min delta = ", min_diff," max delta", max_diff)
    print("arbitrarely pulled 3000th point = ",y1[3000])
    print("signal duration = ",len(theta)/1000)

    
    plt.show()



# Example usage
csv_file   = 'echantillon_test_voie1.csv'  # Replace with your CSV file name or path
csv_file_2 = 'echantillon_test_voie2.csv'
sample_Period = 1 # Period in ms
delta = .012 # distance between lasers in meters
plot_data(csv_file, csv_file_2, sample_Period,delta)

# test = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,11,12,13,14]
# res = average_elements(test,3)
# print(res)