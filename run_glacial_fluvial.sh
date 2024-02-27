#!/bin/bash  

uplift = '0.001 0.0005 0.0001'
initial_grain_size='0.5 0.25 0.1 0.05 0.01'
attrition_rate='0.00001 0.00005 0.0001'
glaciated_switch='True False'

for U in $uplift
do
    for D0 in $initial_grain_size  
    do     
        for a in $attrition_rate
        do
            for glaciation_switch in $glaciated_switch
            do
                python driver.py $U $D0 $a $glaciation_switch
            done
        done
    done
done

