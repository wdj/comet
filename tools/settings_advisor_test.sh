#!/bin/bash
#------------------------------------------------------------------------------

./settings_advisor --help

./settings_advisor --num_node 6000 --platform Titan \
                   --num_vector 28342758 --num_field 882 \
                   --metric_type ccc --num_way 2

./settings_advisor --num_node 6000 --platform Titan \
                   --num_vector 28342758 --num_field 882 \
                   --metric_type czekanowski --num_way 2

./settings_advisor --num_node 3 --platform Titan \
                   --num_vector 450 --num_field 1000 \
                   --metric_type czekanowski --num_way 3 \
                   --sparse no








#------------------------------------------------------------------------------
