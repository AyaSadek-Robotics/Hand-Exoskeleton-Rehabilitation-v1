%% This Function calculates the
%The numerical values for the articualations lengths are considered from 
%Data by ref[7]
  
clc
clear
syms Lpp Lmp Ldp L_tmc L_mcp L_ip th_MCP_abd th_MCP_flex th_PIP th_DIP


DH_table_index=[0 pi/2 0 th_MCP_abd;
              Lpp 0 0 th_MCP_flex;
              Lmp 0 0 th_PIP;
              Ldp  0 0 th_DIP;
];
k=1;
T01=simplify(dhTransform( DH_table_index(k,1),DH_table_index(k,2),DH_table_index(k,3),DH_table_index(k,4) ));k=k+1;
T12=simplify(dhTransform( DH_table_index(k,1),DH_table_index(k,2),DH_table_index(k,3),DH_table_index(k,4) ));k=k+1;
T23=simplify(dhTransform( DH_table_index(k,1),DH_table_index(k,2),DH_table_index(k,3),DH_table_index(k,4) ));k=k+1;
T34=simplify(dhTransform( DH_table_index(k,1),DH_table_index(k,2),DH_table_index(k,3),DH_table_index(k,4) ));k=k+1;



X1=[T04(1,4)];
Y1=[T04(2,4)];
Z1=[T04(3,4)];




[alpha, beta, gamma]=get_orientation(T04(1:3,1:3));





