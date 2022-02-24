##### Energy Efficient VM Placement 


### Sets 

set N;  # set of IP over WDM network nodes -> Core IP router connected to a optical switch
# metro node connected to core node n
set VM;  # set of VM services 
set P; # set of PON networks
set Nm{N} within N;   # set of neighboring of node m
set T;     # should i include to model traffic across different times? Peak traffic and lowest traffic                                                       

##### Questions
# 1 R_c_S ROUTER aggregation ports in the cloud?
# 2 is fog processind node connected to core node s via OLT? 
# 3 Number of router ports used in the fog located in the access network p ?
# 4 Set of neighboring nodes? Multiple sets for nodes in Metro and Access networks?
# 5 Difference between traffic from cloud node s to node d vs traffic from VM v on node s to users in node d

# i<>j for small chars on top

### cloud and fog node parameters
# a fog is consisted of a processing server where a VM can be hosted, along with an ethernet switch, and edge router
# a cloud is said to be hosted on node s if node s is connected to the cloud fog
# bit rate refers to the quantity of data being transferred from one part of the network to the other in a certain time
# bandwidth refers to how much capacity a channel has to send or recieve information

param SW_CB > 0;    #  Cloud switch bit rate.
param SW_CP > 0;     # Cloud switch power consumption.

param SW_MFB > 0;    # Metro fog switch bit rate.
param SW_MFP > 0;    # Metro fog switch power consumption.

param SW_AFB > 0;      # Access fog switch bit rate.
param SW_AFP > 0;      # Access fog switch power consumption.

param SW_R > 0;        # Cloud and fog switch redundancy; Incase a link goes down

param R_CB > 0;  # Cloud router port bit rate.
param R_CP > 0;  # Cloud router port power consumption.

param R_MFB > 0;   # Metro fog router port bit rate.
param R_MFP > 0;   # Metro fog router port power consumption.

param R_AFB > 0;   # Access fog router port bit rate.
param R_AFP > 0;   # Access fog router port power consumption.

param S_P > 0;     # Power consumption of a server.
param S_maxW > 0;    # Maximum workload of a server.

# ratio of the energy used as a whole as compared with the energy used for purely processing
param c > 0;    # Cloud power usage effectiveness. -> ratio of how much energy used for 'real computation' / overhead
param m > 0;    # Metro fog power usage effectiveness.
param a > 0;    # Access fog power usage effectiveness.

# VM parameters
param V >= 0; # number of VMs
param S_v{v in VM} >= 0;   # number of users accessing VM v
param r_v{v in VM} >= 0;  # user download rate of VM v
param L > 0;    # large enough number
param x > 0;    # max number of users served by a single VM replica
param W_v{v in VM} > 0;  # max workload of VM v GHz or as ratio of CPU capacity
param M > 0; # the minimum CPU utilization needed in the absence of load -> when no user is accessing it
param T_v{v in VM} = x * (r[v]);     # traffic resulting from VM replica v, serving the maximum number of users
param W_R_v{v in VM} = (W_v[v] - M) / (T[v]);  # workload per traffic unit, evaluated for VM replica v

### Cloud and fog variables

# hosting cloud
var Cs{s in N} >= 0; # 1 if a cloud is hosted in node s, otherwise 0 -> if node (optical switch and core IP router) host a cloud
var Lsb_C_vs{v in VM, s in N} >= 0;  # v,s = 1 if the cloud hosted in node s hosts a copy of VM v, otherwise v,s = 0 (vm placement in cloud nodes)
var R_C_s{s in N} >= 0;  # Number of router aggregation ports in the cloud hosted in node s.
# number of switches and routers and processing servers
var SW_C_s{s in N} >= 0;   # Number of switches in the cloud hosted in node s.
var S_C_s{s in N} >= 0;  # Number of processing servers in the cloud hosted in node s.

# hosting metro
var F_MF_s{s in N} >= 0;  # Fs = 1 if a fog processing node is hosted in the metro network connected to core node s, otherwise Fs = 0.
var Lsb_MF_vs{v in VM, s in N} >= 0;    # δ(MF) v,s = 1 if the fog processing node hosted inthe metro network connected to node s hosts a
# replica of VM v, otherwise 0
# number of switches and routers in metro fog processing nodes
var R_MF_s{s in N} >= 0;     # Number of router ports used in the fog processing node hosted in the metro network connected to node s
var SW_MF_s{s in N} >= 0;    # Number of switches used in the fog processing node hosted in the metro network connected to node s.
var S_MF_s{s in N} >= 0;    # Number of processing servers in the fog processing node hosted in the metro network connected to node s.

# hosting access
var F_AF_ps{p in P, s in N} >= 0;    # if a fog processing node is built in access network p connected to core node s
var Lsb_AF_vps{v in VM, p in P, s in N} >= 0;   # 1 if the fog processing node in access network p connected to core node s, hosts a replica of VM v
#number of switches, router ports, 
var R_AF_vps{v in VM, p in P, s in N} >= 0;     # Number of router ports used in the fog processing node located in the access network p connected to core node s
var R_AF_ps{p in P, s in N} >= 0;    # Number of router ports used in the fog processing node located in the access network p connected to core node s.
var SW_AF_ps{p in P, s in N} >= 0;    # Number of switches used in the fog processing node located in access network p connected to core node s.
var S_AF_ps{p in P, s in S} >= 0;     # Number of processing servers in the fog processing node located in the access network p connected to core node s.

# VM variables
var W_CR_vs{v in VM, s in N} >= 0;   # workload of VM replica v hosted in node S
var W_C_s{s in N} >= 0; # total workload of cloud hosted in node s
var D_C_vsd{v in VM, s in N, d in N} >= 0; # traffic flow from VM replica v to, hosted in a cloud node s, to users in node d
var L_sd{s in N, d in N} >= 0;      # should I do for neighboring nodes ie. {s in N, d in Nm[s]}? 
var W_MFR_vs{v in VM, s in N} >= 0;   # workload of the VM replica v hosted in fog processing node s 
var W_MF_s{s in N} >= 0;  # total workload of the metro fog processing node located in core node S
# does this mean that metro processing node is in a core node 
var D_MF_vs{v in VM, s in N} >= 0;  # traffic from VM replica v hosted connected to node S
var W_AFR_vps{v in VM, p in P, s in N} >= 0;  # workload of VM replica v hosted in fog processing node located in access network p connected to core node S
var W_AF_ps{p in P, s in N} >= 0;  # total workload of the fog processing node located in the access netwokr p connected to core node S
var D_AF_vps{v in VM, p in P, s in N} >= 0; #  traffic flow from the VM replica v hosted in the fog processing node located in the access netwokr p connected to core node S


# Access network parameters
param A_p{p in P} >= 0; # average broadband data rate in PON p 
param O_v{v in VM} >= 0;   # ratio of traffic due to VM c to the total PON traffic # VM v traffic / 
# total traffic - all traffic of users (emaiil, video) 
param OLT_B_pd{p in P, d in N} >= 0;

param U_vpd{v in VM, p in P, d in VM} = (OLT_B_pd[p,d] / A_p[p]) * O_v[v];
# total OLT capacity / average data rate (phi = how many users are connected to PON)= number of users
# total number of users it can support


param OLT_N_pd{p in P, d in N} >= 0;
param OLT_P >= 0;
param D_vpd{v in VM, p in P, d in N} = U_vpd[v, p, d]; # 
param ONU_N_pd{p in P, d in N} >= 0;
param ONU_P > 0;  # power consumption of an ONU 
parma n > 0; # network power usage effectiveness 

# Metro network parameters
param R_MB > 0;  # metro router bit rate
param R_MP > 0;  # metro router power consumption
param R_MR > 0;   # metro router redundancy
param SW_MB > 0;  # metro ethernet switch bit rate
param SW_MP > 0;   # metro ethernet power consumption

# Metro network variables 
var R_M_s{s in N} >= 0;   # number of router ports in metro netwokr connected to node S
var SW_M_s{s in N} >= 0;     # number of ethernet switches in metro network connected to node s 

# core network parameters
param R_p > 0; # core router port power consumption 
param t_p > 0; # transponder power consumption
param e_P > 0; # EDFA power consumption
param SW_P_s{s in N} > 0;  # optical power consumption in node S
param G_P > 0; # regenerator power consumption
param W > 0; # number of wavelenghts per fibre
param W_B > 0; # wavelengh data rate 
param S > 0; # maximum span distance between two EDFAs in km 
param D_mn{m in N, n in N} >= 0; # distance in km between node pair (m,n)
param A_mn{m in N, n in Nm[m]}= ((D_mn[m,n]/S) - 1) +2 ;   # No.Of EDFAs bet WDM nodes, 2 is for post and pre amplification and the -1 is to ensure we dont add EDFA for first 80km
param G_mn{m in N, n in Nm[m]} = ((D_mn[m,n]/R) - 1);

# core network variables
var C_ij{i in N,j in N,t in T:i<>j}>=0;                                 # No. of WL in virtual link i-j																									
var W_mn {m in N, n in Nm[m], t in T}>=0;                                # total No. of WLs in a physical link																									
var R_AC_s{s in N} >= 0; # number of router ports in node s that aggregate traffic from/to clouds
var R_AE_d{d in N} >= 0;  # number of router ports in node d that aggregate the traffic from/to metro routers
var F_mn{m in N, n in Nm[m]} >= 0;   # number of fibres on a physical link (m,n)
var Lsd_ij {s in N,d in N,i in N,j in N,t in T:s<>d}>=0;               # portion of s-d demand passing thru i-j virtual link																									
var Wij_mn {i in N,j in N, m in N, n in Nm[m], t in T:i<>j}>=0;         # portion of WLs in i-j passing thru m-n																									



# Objective -> minimize the total power consumption in cloud to fog architecture
# Core + Metro + PON + CLOUD + MF + AF 
# CLOUD 
minimize CF :

# power consumption of the cloud servers
c * sum{s in N}(S_C_s[s]*S_P) +  # PUE * iterate through all nodes -> multiply processing server by power consumption
# power consumption of cloud switches and routers
c * (sum{s in N}((SW_C_s[s]*SW_R*SW_CP) + (R_C_s[s]+R_CP)) +

# Metro Fogs
# power consumption of the metro servers
m * sum{s in N}((S_MF_s[s]*S_P)) +
# power consumption of the swithces and routers
m * sum{s in N}((SW_MF_s[s]*SW_R*SW_MFP) + (R_MF_s[s] * R_MFP)) + 

# Access Fog
# power consumption of the access servers
a * sum{s in N}(sum{p in P} S_AF_ps[p,s] * S_P) + 
# power consumption of switches and routers
a * (sum{s in N}(sum{p in P}((SW_AF_ps[p, s] * SW_R + SW_AFP) + R_AF_ps[p,s]) + 

# PON 
# contains the ONU -> end user interface, OLT -> connects multiple ONUs
# power consumption of an OLT
n * (sum{p in P}(sum{d in N}(OLT_P * OLT_N_pd[p, d]))) + 

# power consumption of ONUs
n * (sum{p in P}(sum{d in N}(ONU_P * ONU_N_pd[p, d]))) + 

# Metro 
n * (sum{s in N}R_M_s[s] * R_MR * R_MP) + 
n * (sum{s in N}SW_M_s[s] * SW_MP) + 

# Core network
n * (sum{s in N}(R_P * R_AC_s[s]) 
+ sum{d in N}(R_P * R_AE_d[d]) 
+ sum{m in N}(sum{n in Nm[m]}R_P * W_m_n[m,n])) + 

n * (sum{m in N}(sum{n in Nm[m]}) t_p * W[m,n] + 

n * sum{m in N}(sum{n in Nm[m]}(e_P * F_mn[m,n] * A_mn[m, n]) + 

n * (sum{s in N}SW_P_s[s]) + 

n * (sum{m in N}(sum{n in Nm[m]:m<>n})G_P * G_mn[m, n] * W_mn[m, n]);
#### Constraints


# Placing VM in cloud constraints
# serving VM demand constraints
subject to Cons17{v in VM}:  
# traffic flow from vm v to users in PON network p connected to core node d
# sum of p in P and sum d in N
# each core node -> we consider all the 
# traffic between users in PON p for VM v over all the destination
sum{p in P}(sum{d in N})D_vpd[v,p,d] = (sum{s in N}(sum{d in N}D_C_vsd[v,s,d])) + (sum{s in N}D_MF_vs[v,s]) + 
(sum{p in P}(sum{s in N}D_AF_vps[v,p,s]));
# ensure that the users are served by VM v
# traffic from VM v -> downlink traffic, uplink traffic (small fraction of downlink)
# traffic from two ways 
# placed in access fogs metro fogs core nodes 

subject to Cons18{s in N, v in VM}:
L * (sum{d in N}(D_C_vsd[v,s,d])) >= Lsb_C_vs[v,s];

subject to Cons19{s in N, v in VM}:
(sum{d in N}(D_C_vsd[v,s,d])) >= L * Lsb_C_vs[v,s];

# placing VM in metro fog constraint
subject to Cons20 {s in N, v in VM}:
D_MF_vs[v,s] >= Lsb_MF_vs[v,s];

subject to Cons21 {s in N, v in VM}:
D_MF_vs[v,s] <= L * Lsb_MF_vs[v, s];

# placing VM in access fog constraints 
subject to Cons22{s in N, v in VM, p in P}:
D_AF_vps[v, p, s] >= Lsb_AF_vps[v, p, s];

subject to Cons23{s in N, v in VM, p in P}:
sum{v in VM D_AF_vps[v, p, s] >= L * Lsb_AF_vps[v, p, s];

# Cloud location constraint
subject to Cons24{s in N}:
sum{v in VM} Lsb_C_vs[v, s] >= C_s[s];

subject to Cons25{s in N}:
sum{v in VM} Lsb_C_vs[v, s] >= L * C_s[s];

# Metro fogs location constraints 
subject to Cons26{s in N}:
sum{v in VM} Lsb_MF_vs[v,s] >= F_MF_s[s];

subject to Cons27{s in N}:
sum{v in VM} Lsb_MF_vs[v,s] >= L * F_MF_s[s];

# Access fog location constraints
subject to Cons28{s in N}:
sum{v in VM} Lsb_AF_vps[v, p, s] >= F_AF_ps[p, s];

subject to Cons29{s in N}:
sum{v in VM} Lsb_AF_vps[v, p, s] >= L * F_AF_ps[p, s];

# Cloud and fog workload constraints
# VM replica workload in a cloud fog
subject to Cons30{v in VM, s in N}: 
W_CR_vs[v, s] = Lsb_C_vs[v, s] * W_v[v];

subject to Cons31{v in VM, s in N}: 
W_CR_vs[v, s] = (((sum{d in N}D_C_vsd[v, s, d]) / (r_v[v] * x)) * M * Lsb_C_vs[v,s]) + 
(W_R_v[v] * (sum{d in N}D_C_vsd[v, s, d]); 

subject to Cons32{s in N}:
W_C_s[s] = sum{v in VM} W_CR_vs[v,s];

subject to Cons33{v in VM, s in N}:
W_MFR_vs[v,s] = Lsb_MF_vs[v,s] * W_v[v];

subject to Cons34{v in VM, s in N}:
W_MFR_vs[v,s] = (((D_MF_vs[v,s])/(r_v[v] * x)) * M * Lsb_MF_vs[v,s]) + (W_R_v[v] * D_MF_vs[v,s]);

subject to Cons35{s in N}: 
W_MF_s[s] = sum{v in VM} W_MFR_vs[v,s];

subject to Cons36{v in VM, s in N, p in P}: 
W_AFR_vps[v,p,s] = Lsb_AF_vps[v,p,s] * W_v[v];

subject to Cons37{v in VM, s in N, p in P}:
W_AFR_vps[v,p,s] = (((S_AF_vps[v,p,s])/(r_v[v]*x)) * M * Lsb_AF_vps[v,p,s]);

subject to Cons38{s in N}: 
W_AF_ps[p,s] = sum{v in VM} W_AFR_vps[v,p,s];

# Number of servers in cloud and fog constraints:

subject to Cons39{s in N}: 
S_C_s[s] >= W_C_s[s] / S_maxW;

subject to Cons40{s in N}:
S_MF_s[s] >= W_MFR_vs[v, s] / S_maxW;

subject to Cons41{s in N, p in P}: 
S_AF_ps[p,s] >= W_AF_ps[p,s] / S_maxW;

# Number of router ports and switches in cloud and fog
subject to Cons42{s in N}:
R_C_s[s] >= ((sum{s in VM}(sum{d in N}D_C_vsd[v,s,d])) / R_CB;

subject to Cons43{s in N}:
SW_C_s[s] >= ((sum{v in VM}(sum{d in N})D_C_vsd[v,s,d])) / SW_CB;

subject to Cons44{s in N}: 
R_MF_s[s] >= ((sum{v in VM}D_MF_vs[v,s]) / R_MFB);

subject to Cons45{s in N}: 
SW_MF_s[s] >= ((sum{v in VM}D_MF_vs[v,s]) / SW_MFB);

subject to Cons46{s in N, p in P}: 
R_AF_ps[p,s] >= (sum{v in VM}D_AF_vps[v, p, s] / R_AFB);

subject to Cons47{s in N, p in P}: 
SW_AF_ps[p, s] >= (sum{v in VM}D_AF_vps[v,p,s]/SW_AFB);

# Number of metro router ports and ethernet switches in metro network constraints
subject to Cons48{s in N}:
R_M_s[s] >= (sum{v in VM}(sum{s in N} D_C_vsd[v, s,d] + sum{v in VM}D_MF_vs[v,s])) / R_MB;

subject to Cons49{s in N}:
SW_M_s[s] >= (sum{v in VM}(sum{s in N} D_C_vsd[v, s,d] + sum{v in VM}D_MF_vs[v,s])) / SW_MB;

subject to Cons50{s,d in N}: 
sum{v in VM} D_C_vsd[v,s,d];

##########flow conservation in IP layer...Ref (2)
s.t. Cons2{s in N, d in N, i in N:s<>d}:
sum{j in N:i<>j}Lsd_ij[s,d,i,j] - sum{j in N:i<>j}Lsd_ij[s,d,j,i]=
if i=s then Lsd[s,d] else if i=d then -Lsd[s,d] else 0;

#####accumilated demands portions on a virtual link does not exceed its capacity..(3)
s.t. Cons3{i in N, j in N:i<>j}:
sum{s in N}(sum {d in N:s<>d} Lsd_ij[s,d,i,j])<=Cij[i,j]*B;

##### total No. of ports does not exceed the limit..(4)
s.t. Cons4{i in N}:sum{j in N:i<>j}Cij[i,j]+Aggrp[i]<=rpMax[i];

#####total No. of ports does not exceed the limit..(5)
s.t. Cons5{j in N}:sum{i in N:i<>j}Cij[i,j]+Aggrp[j]<=rpMax[j];

#####flow conservation in WDM layer...(6)
s.t. Cons6{i in N, j in N, m in N:i<>j}:sum{n in Nm[m]}Wij_mn[i,j,m,n]-sum{n in Nm[m]}Wij_mn[i,j,n,m]=
if m=i then Cij[i,j] else if m=j then -Cij[i,j] else 0;

#####the accumilated portions of virtual links wavelenghs on a physical link does not exceed its capacity..(7)
s.t. Cons7{m in N, n in Nm[m]}:
sum{i in N}(sum{j in N:i<>j} Wij_mn[i,j,m,n])<=W*Fmn[m,n];

## constrains for Wmn so they dont be zero..
s.t. Cons8{m in N, n in Nm[m]}:
sum{i in N}(sum{j in N:i<>j} Wij_mn[i,j,m,n])=Wmn[m,n];

## constrains of calculations of used data aggregation ports
s.t. Cons9{s in N}: Aggrp[s] =sum{d in N:d<>s}Lsd[s,d]/B;


