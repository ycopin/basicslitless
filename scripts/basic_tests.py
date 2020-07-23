#### Pour le rapport :



from starsanalysis import *

# config = configuration() #default configuration initialisation
# config.set_wcs(0,0) #wcs centered on 0,0
# 
# star = star_class(0,0,5,config)
# star.set_orders([-2,-1,0,1,2],config)
# 
# fig,ax = plt.subplots()
# for n in star.order:
#     ax.add_patch(PolygonPatch(star.order[n], color = np.random.choice(range(256), size=3)/256, alpha = 0.7))
# 
# ax.set_xbound(0,2048)
# ax.set_ylim(900,1100)
# ax.axis('equal')
# plt.show()


config = configuration("C:/Users/Sadek/basicslitless/config/auxtel.ini")

cstar = star_class(0,0,5,config)

maglim = 5 + 5

cstar.set_orders([0,1], config, maglim, 1)


cstar.translate(-cstar.x, -cstar.y)
L = []
Mag = [1, 2.5,4 ,5]
n = len(Mag)
for i in Mag:
    L.append(star_class(0,0,5+i,config))
    L[-1].set_orders([0,1],config,maglim,1)
    L[-1].translate(-L[-1].x,-L[-1].y)
    
fig, ax = plt.subplots(1,n)

for i in range(n):
    ax[i].add_patch(PolygonPatch(cstar.all_orders, color = 'red', alpha = 0.5))
    ax[i].add_patch(PolygonPatch(L[i].all_orders, color = 'blue'))
    ax[i].set_xlim(-100, 2000)
    ax[i].set_ylim(-20,20)
    ax[i].set_xlabel('Position X (px)')
    if i == 0:
        ax[i].set_ylabel('Position Y (px)')
    ax[i].set_title('Superposition : {:.2f} %'.format(total_overlap(cstar,L[i],0)))
    
plt.show()