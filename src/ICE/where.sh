#!/bin/bash

#grep -ri 'v_i' *90 | grep -Ev 'v_il|v_ip|sv_i|v_ice|bv_i|pv_i|v_i_old|v_i_t|v_i_t_1d|v_i_t_2d|zv_i|z1_v_i|dv_i|v_ial'

grep -ri 'v_s' *90 | grep -Ev 'v_sl|v_sp|sv_s|v_sce|bv_s|pv_s|v_s_old|v_s_t|v_s_t_1d|v_s_t_2d|zv_s|z1_v_s|dv_s|v_sal|v_str|v_sub'
