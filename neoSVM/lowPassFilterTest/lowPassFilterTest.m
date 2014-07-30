% lowPassFilter Test
function lpData = lowPassFilterTest(data,sr)

filter_n=256; % 100*sr
lowP=16; %cutoff frequency
highP=3;
lowpass=fir1(filter_n,[highP*2/sr lowP*2/sr]);
lpData=filter(lowpass,1,data);