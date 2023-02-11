import numpy as np
import matplotlib.pyplot as plt

def periodic(coord, box):
        if abs(coord) > 0.5 * box: 
            return coord - np.sign(coord) * box
        return coord 
    
class Box():
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    def periodic_correct(self, xb, yb, zb):
        xb = periodic(xb, self.x)
        yb = periodic(yb, self.y)
        zb = periodic(zb, self.z)
        return xb, yb, zb
    
def rnd_vector(length_bond=1.0):
    v = np.random.uniform(-1.0, 1.0, 3)
    v /= np.sqrt(np.sum(v**2)) * length_bond 
    return v

def chain(na, nb, box, period, length_side_chain): # period after which side chains appear
    list_coord = list()
    list_bond = list()
    v = np.random.uniform(-1.0, 1.0, 3)
    x = v[0] * box.x/2
    y = v[1] * box.y/2
    z = v[2] * box.z/2
    list_coord.append([x, y, z, 1])
    for i in range(1,na):
        v = np.random.uniform(-1.0, 1.0, 3)
        x += v[0]
        y += v[1]
        z += v[2]
        x, y, z = box.periodic_correct(x, y, z)
        list_coord.append([x, y, z, 1])
        if len(list_coord) % period == 0: #initiating the addition of beads to the side chain na, as a continuation of the main chain
            v = np.random.uniform(-1.0, 1.0, 3)
            x_side = v[0] + x
            y_side = v[1] + y
            z_side = v[2] + z
            x_side, y_side, z_side = box.periodic_correct(x_side, y_side, z_side)
            list_coord.append([x_side, y_side, z_side, 1])
            for _ in range(1,length_side_chain):
                v = np.random.uniform(-1.0, 1.0, 3)
                x_side += v[0]
                y_side += v[1]
                z_side += v[2]
                x_side, y_side, z_side = box.periodic_correct(x_side, y_side, z_side)
                list_coord.append([x_side, y_side, z_side, 1])  #after that return to chain frame
    
    for i in range(1,nb+1):  
        v = np.random.uniform(-1.0, 1.0, 3)
        x += v[0]
        y += v[1]
        z += v[2]
        x, y, z = box.periodic_correct(x, y, z)
        list_coord.append([x, y, z, 2])   
        if len(list_coord) % period == 0:#initiating the addition of beads to the side chain nb
            v = np.random.uniform(-1.0, 1.0, 3)
            x_side = v[0] + x
            y_side = v[1] + y
            z_side = v[2] + z
            x_side, y_side, z_side = box.periodic_correct(x_side, y_side, z_side)
            list_coord.append([x_side, y_side, z_side, 2])
            for _ in range(1,length_side_chain):
                v = np.random.uniform(-1.0, 1.0, 3)
                x_side += v[0]
                y_side += v[1]
                z_side += v[2]
                x_side, y_side, z_side = box.periodic_correct(x_side, y_side, z_side)
                list_coord.append([x_side, y_side, z_side, 2])
    return list_coord      

if __name__ == '__main__':
    box = Box(15, 15, 15)
    na = 10
    nb = 4
    period = 3
    length_side_chain = 3
    n_chain = 50
    n_solent = 15 ** 3 * 3 - n_chain * (na + nb)
    n = (na + nb) * n_chain + n_solent
    n_bonds = (na + nb - 1) * n_chain
    fcoord = open('COORD', 'w')
    fbonds = open('BONDS', 'w')
    fcoord.write(f'num_atoms {n} box_size {box.x} {box.y} {box.z}\n')
    fbonds.write(f'num_bonds {n_bonds} num_atoms {n} box_size {box.x} {box.y} {box.z}\n')
    temp = 0
    for _ in range(n_chain):
        coord = chain(na=na, nb=nb, box=box, period=period, length_side_chain=length_side_chain) #create a chain
        for i, w in enumerate(coord):
            fcoord.write(f'{i+1+temp: <10} {w[0]: <25} {w[1]: <25} {w[2]: <25} {w[3]}\n'.format(*w)) #write the coordinates to the file
        length_na = na + (na // period) * length_side_chain
        for i in range(1, length_na):   #writning bonds for na
            if i % (period + length_side_chain) == 0: # if 
                fbonds.write(f'{i + temp - length_side_chain} {i + 1 + temp}\n')
                continue
            fbonds.write(f'{i + temp} {i + 1 + temp}\n')
        fbonds.write(f'{length_na + temp} {length_na + 1 + temp}\n')
        for i in range(1,nb + (nb // period) * length_side_chain):
            if i % (period + length_side_chain) == 0:
                fbonds.write(f'{i + temp - length_side_chain + length_na} {i + 1 + temp + length_na}\n')
                continue
            fbonds.write(f'{i + temp + length_na} {i + 1 + temp + length_na}\n')    
        temp += na + nb + (na // period) * length_side_chain + (nb // period) * length_side_chain
    fbonds.close()  
    x_s = np.random.uniform(-box.x / 2, box.x / 2, n_solent)
    y_s = np.random.uniform(-box.y / 2, box.y / 2, n_solent)
    z_s = np.random.uniform(-box.z / 2, box.z / 2, n_solent)    
    for x,y,z in zip(x_s, y_s, z_s):
        temp += 1
        fcoord.write(f'{temp: <10} {x: <25} {y: <25} {z: <25} 3\n'.format(x, y, z, temp))
    fcoord.close()

    # fig, ax = plt.figure(), plt.axes(projection ='3d')
    # diblock = chain(na=na, nb=nb, box=box, period=period, length_side_chain=length_side_chain)
    # # print(diblock)
    # x, y ,z ,c = [], [], [], []
    # for bead in diblock:
    #     print(bead)
    #     x.append(bead[0])
    #     y.append(bead[1])
    #     z.append(bead[2])
    #     c.append(bead[3]) 
    # ax.scatter(x, y, z, "-o", c = c)
    # ax.plot3D(x, y, z, "-")
    # plt.show()