function [kfx,kfy,kfz] = CROSS2(di_x,di_y,di_z,dj_x,dj_y,dj_z)

kfx =  di_y.* dj_z - di_z.* dj_y;
kfy =-(di_x.* dj_z - di_z.* dj_x);
kfz =  di_x.* dj_y - di_y.* dj_x;