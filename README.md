# LTE_EVA_ChannelwithMIMO
## Basic Simulation Conditions
- Band: 20MHz
- Carrier frequence: 2GHz
- Channel: Extended Vehicular A Model(EVA channel)
- Speed: 60km/h, 250km/h
- Total Data Source Volume: >1e7

## Design Block Diagram
### Transmitter System
- Signal generated
- Convolution/Turbo encoding
- Random interleaver
- Modulation(BPSK, QPSK, 16QAM)
- Signal divided
- Add pilots
- Serial-to-parallel conversion
- IDFT
- Parallel-to-serial conversion
- add CPs
- MIMO(V-blast/STBC)

### channel
- Rayleigh channel
- AWGN channel

### Receiver System
- Remove CPs
- Serial-to-parallel conversion
- DFT
- Parallel-to-serial conversion
- Channel estimate(LS, MMSE algorithm)
- MIMO receiver(ML algorithm)
- Time synchronization
- Demodulation(BPSK, QPSK, 16QAM)
- Deinterleaver
- Convolution/Turbo decoding(Viterbi Algorithm)
- Signal received

Here, I can provide a figure with you using **Simulink and Model-Based Design**, but sadly, it's just a **_SISO_** design diagram but **_NOT MIMO_**. However, I do not think they have differences in nature though MIMO is somewhat more complicated.
![SISO](https://user-images.githubusercontent.com/40145471/129459558-5a2235ca-f3e4-4fc5-bd2f-f4d8042dacbd.png)
Fig. A SISO System designed by Simulink

## Conclusion
### Parameters
After carefully and prudentially comparing, the values of parameters are as follows:
| Parameter | Value |
|:-----------:|:-------:|
|Interleaver Depth| 846 |
|Source Encoding| Turbo|
|CP Length| 20 |
|Pilot Gap|20|

I am sorry for not writing down why I am choosing these parameters here since it was really complicated and I would waste plenty of time explaining it, but if you want to know, **feel free to contact me at _qshan.yuezhao@gmail.com_!** or you can explore [LTE_EVA_ChannelwithMIMO/LTE_test](https://github.com/DaGuanYuan/LTE_EVA_ChannelwithMIMO/new/master/LTE_test) by yourself (it's really a tough work :pensive:).

### Results
By Comprehesive consideration of **Communication Quality(Average BER), Transmission Rate, and Extra cost(produced by CP, pilot, encoding....)**, here preferred modulation methods under different communication circumstances are recommended and reference transmission speed is given as below:

|SNR|Preferred Modulation Methods|Transmission Rate|
|:-:|:--------------------------:|:---------------:|
|<=10dB|BPSK|7.32Mbit/s|
|10~20dB|QPSK|14.64Mbit/s|
|>=30dB|16QAM|29.28Mbit/s|
