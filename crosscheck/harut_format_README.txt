info on data file columns from emails from Harut

> I put a text file in
> /work/clas12/avakian/mc/mcjune2019/T-1.00_S-1.0/clasdis/cooked6b20/dihad1000.dis.0000.nrad.dat.evio.hipo.txt
>
> with the following structure (all generated quantities) for events
> with 2 pions
>
> event#,pip_energy,pip_theta,pip_phi,pip_fstar,pip_pt,pip_zi
> event#,pim_energy,pim_theta,pim_phi,pim_fstar,pim_pt,pim_zi
> event#,P_h-phistar,P_h-P_T,phi_Rt1,phi_Rt2
>
> R_T1 calculated from R_T=0.5(Pt1-Pt2)
> R_T2 calculated from R_T2=z2*Pt1-z1*Pt2


UPDATE:
It has the event number from the RUN::config bank (header.json), and in
addition the list of variables (all generated) extended

event#,pip_energy,pip_theta,pip_phi,pip_fstar,pip_pt,pip_zi,pip_mis_max,pip_xf,pip_et,pip_etbr
(last 2 rapitity in the CM and Breit frames)

definitions:

pip_phi -> Lab phi of pi+

pip_theta -> Lab theta of the pi+

pip_fstar -> pi+ \phi  (between pi+ gamma* and e\gamma* planes (standard
Trento convention)

pip_mis_max ->missing mass of the e'\pi+X system

pip_et,pip_etbr --> pseudorapidity in gamma*-N CM and Breit Frame.
(formulas were checked with http://arxiv.org/pdf/1611.10329.pdf)

P_h-phistar --> this is the phi_star of the dihadron (sum of two
hadrons) defined the same way as for a single hadron
  
UPDATE:
Now the third row for each event has in addition the theta and also phi angle between VM
(dihadron)  and pi+.
