#!/bin/bash
# Run the 'cobold.job' with a nohup background and then view with htop
rm rhd.stop
nohup cobold.job&
htop
