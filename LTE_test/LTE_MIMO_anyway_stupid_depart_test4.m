%LTE_7:交织器和译码器的选择研究
%% 初始条件和信号生成%%%%%%%%%%%%%%%%%%%%%%%%%%
num_bit = 1.28e4;
Nsamp = 4;
M = 4;
k = log2(M);
Rs = 20*(10^6);    %信道带宽
Rb = Rs*k;    %bit rate 
Tb = 1/Rb;   
Ts = 1/Rs;
rng('default');
rng(67);    %seed = 67
msg = randi([0,1],1,num_bit);  %generate messege
SNR_awgn  = 20;%awgn snr
%%MIMO%%
numTx = 2;
numRx = 2;
chanSRate = Rs;
% PathDelays = [0,30,150,310,370,710,1090,1730,2510]*1e-9 ;
% PathGains =  [0,-1.5,-1.4,-3.6,-0.6,-9.1,-7.0,-12.0,-16.9];
PathDelays = [0]*1e-9 ;
PathGains =  [0];
velocity = 60/3.6;      %velocity
f = 2*10^9 ;     %载频
c = 3*10^8;     %light speed
Doppler = velocity*f/c;     %Maxinum Doppler shift
% Doppler = 100;



%% %%%%%%%卷积编码%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
constlen = 7;   %constraint length
codegen = [171, 173];   % polynomials generator
tblen = 5*constlen;     %traceback length
trellis = poly2trellis(constlen, codegen);
dspec = distspec(trellis, 7);
msg_conv = convenc(msg, trellis);   %convolutional encode
msg_turbo_cascade = msg_conv;%Here msg_conv_re could be both turbo encode and convolution encode
G = length(msg_conv)/length(msg);    %code generator gain of the convolutional encode


%% Turbo编码 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% frmlen = num_bit;     %The frame length ->may not be useful
% intrlvrIndices = randperm(num_bit);
% Turbo_encoder = comm.TurboEncoder(...
%     'TrellisStructure', poly2trellis(4, [13 15 17], 13),...
%     'InterleaverIndices', intrlvrIndices);  %This is the Turbo encoder
% msg_turbo = step(Turbo_encoder, transpose(msg));
% msg_turbo_cascade = transpose(msg_turbo);
% msg_conv = msg_turbo_cascade;  %Here msg_conv_re could be both turbo encode and convolution encode
% G = length(msg_turbo_cascade)/length(msg);  %code generator gain of the turbo encode

% msg_conv = msg;     %test
%% 行列交织 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% row_intrlv = 1;     %interleave depth
% column_intrlv = 128;   %interleave width
% group_num = fix(length(msg_conv)/(row_intrlv*column_intrlv));  %row_intrlv*column_intrlv bits as a group;group number
% more_num1 = length(msg_conv) - row_intrlv*column_intrlv*group_num;
% rest_num1 = row_intrlv*column_intrlv* (group_num + 1) - length(msg_conv);
% group_num = group_num + 1;
% msg_conv(1, group_num*row_intrlv*column_intrlv) = 0;    %补0
% msg_conv_re = transpose(reshape(msg_conv,row_intrlv*column_intrlv,group_num));     %group by group
% msg_intrlv = zeros(group_num,row_intrlv*column_intrlv);     %initial->distribute caches
% for i = 1:group_num
%     msg_intrlv(i,:) = matintrlv(msg_conv_re(i,:), row_intrlv, column_intrlv);    %100*50 matrix interleave
% end
% msg_intrlv_cascade = transpose(reshape(transpose(msg_intrlv), 1, size(msg_intrlv,1)*size(msg_intrlv,2)));   %interleave completed!!!!!!!!!!    

% msg_intrlv_cascade = msg'; %test

%% 随机交织%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
intrl_frm = 1000;
intrl_rand_ts =fix(length(msg_conv) / intrl_frm);   %全都补零部分情况了！！！！
more_num1 = length(msg_conv) - intrl_frm*intrl_rand_ts;
rest_num1 = intrl_frm*(intrl_rand_ts + 1) - length(msg_conv);  %补的0的个数
intrl_rand_ts = intrl_rand_ts + 1;
msg_turbo_cascade(1, intrl_frm*(intrl_rand_ts)) = 0;  %无论如何都补0，反正都去掉
msg_conv_re = transpose(reshape(msg_turbo_cascade, intrl_frm, intrl_rand_ts));%here msg_conv_re could be both turbo encode and convolution encode
msg_intrlv = zeros(intrl_rand_ts , intrl_frm);
for i = 1:intrl_rand_ts
    msg_intrlv(i,:) = randintrlv(msg_conv_re(i,:), i) ;     %i is the random seed, 我真是个天才=-=这样每一个的随机过程都不一样
end
msg_intrlv_cascade = transpose(reshape(transpose(msg_intrlv), 1, size(msg_intrlv,1)*size(msg_intrlv,2)));  %interleave completed!!!!!!!!!!!



%% QPSK调制%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg_syb = transpose(reshape(msg_intrlv_cascade, size(msg_intrlv_cascade, 2)*k, size(msg_intrlv_cascade, 1)/k));     %IQ, get symbols
msg_de = bi2de(msg_syb);    %decimal symbols
grayencode = bitxor(0:M-1, floor((0:M-1)/2));   %generate gray list by order,[0,1,3,2]
msg_gry_enc = grayencode(msg_de+1);     %gray encode:0->0(The first),1->1(The second)，2->3(The third),3->2(The forth)，so by order,msg_enc+1
msg_qpsk = pskmod(msg_gry_enc,M, pi/4);    %QPSK modulation

%% 分成两路%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg_qpsk_div = reshape(msg_qpsk, length(msg_qpsk)/2,2);     %把信号分给两个天线
msg_qpsk1 = transpose(msg_qpsk_div(:,1));      %天   线1
msg_qpsk2 = transpose(msg_qpsk_div(:,2));      %天线2
 
%% 串并变换%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_carrier = 128;     %子载波数量
len_cp = 12;    %length of cp

%% 为天线1添加导频1 天线2添加导频0—>估计h11,h21%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pilotfrequency = 3;     %导频间隔
pilot_num = fix(length(msg_qpsk)/(2*pilotfrequency));    %导频数量
rest_num2 = pilotfrequency*(pilot_num + 1) - length(msg_qpsk)/2;  %单个导频补零
pilot_num =pilot_num + 1;
msg_qpsk1(1, pilot_num * pilotfrequency) = 0;    %补0
msg_qpsk2(1, pilot_num * pilotfrequency) = 0;    %补0

msg_pilot1 = transpose(reshape(msg_qpsk1, pilotfrequency,length(msg_qpsk1)/pilotfrequency));   %trick, 添加导频方便而已
msg_pilot2 = transpose(reshape(msg_qpsk2, pilotfrequency,length(msg_qpsk2)/pilotfrequency));   %trick, 添加导频方便而已

msg_pilot_aft1 = zeros(length(msg_qpsk1)/pilotfrequency,pilotfrequency+1);    %initialize
msg_pilot_aft2 = zeros(length(msg_qpsk2)/pilotfrequency,pilotfrequency+1);    %initialize

for i = 1:length(msg_qpsk1)/pilotfrequency
   msg_pilot_aft1(i,:) = [1, msg_pilot1(i,:)] ;   %天线1，插入导频1
   msg_pilot_aft2(i,:) = [0, msg_pilot2(i,:)] ;   %天线2，插入导频0
end

msg_pilot_aft_casd10_1 = reshape(transpose(msg_pilot_aft1), 1, size(msg_pilot_aft1,1)*size(msg_pilot_aft1,2));     %转成行
msg_pilot_aft_casd10_2 = reshape(transpose(msg_pilot_aft2), 1, size(msg_pilot_aft2,1)*size(msg_pilot_aft2,2));     %转成行
% ab_c : a:第一个天线导频，b:第二个天线导频，c:第c个天线
% msg_pilot_aft_casd10 = cat(2, msg_pilot_aft_casd1,msg_pilot_aft_casd2);
%% 为天线1添加导频0 天线2添加导频1—>估计h12,h22%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

msg_qpsk1(1, pilot_num * pilotfrequency) = 0;    %补0
msg_qpsk2(1, pilot_num * pilotfrequency) = 0;    %补0

msg_pilot1 = transpose(reshape(msg_qpsk1, pilotfrequency,length(msg_qpsk1)/pilotfrequency));   %trick, 添加导频方便而已
msg_pilot2 = transpose(reshape(msg_qpsk2, pilotfrequency,length(msg_qpsk2)/pilotfrequency));   %trick, 添加导频方便而已

msg_pilot_aft1 = zeros(length(msg_qpsk1)/pilotfrequency,pilotfrequency+1);    %initialize
msg_pilot_aft2 = zeros(length(msg_qpsk2)/pilotfrequency,pilotfrequency+1);    %initialize

for i = 1:length(msg_qpsk1)/pilotfrequency
   msg_pilot_aft1(i,:) = [0, msg_pilot1(i,:)] ;   %天线1，插入导频1
   msg_pilot_aft2(i,:) = [1, msg_pilot2(i,:)] ;   %天线2，插入导频0
end

msg_pilot_aft_casd01_1 = reshape(transpose(msg_pilot_aft1), 1, size(msg_pilot_aft1,1)*size(msg_pilot_aft1,2));     %转成行
msg_pilot_aft_casd01_2 = reshape(transpose(msg_pilot_aft2), 1, size(msg_pilot_aft2,1)*size(msg_pilot_aft2,2));     %转成行
% msg_pilot_aft_casd01 = cat(2, msg_pilot_aft_casd1,msg_pilot_aft_casd2);     %用01传数据

%% 串并变换%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group_prl = fix(length(msg_pilot_aft_casd01_1)/num_carrier);   %number of the groups which are divided in serial-to-parallel conversion
more_num3 = length(msg_pilot_aft_casd01_1) - num_carrier*group_prl; %the same as interleaving
rest_num3 = num_carrier*(group_prl + 1) - length(msg_pilot_aft_casd01_1);  %第二次补0个数
rest_num = rest_num2 + rest_num3 ;   %总补零的个数
group_prl = group_prl + 1;

msg_pilot_aft_casd01_1(1, num_carrier * group_prl) = 1;   %补0是为了把不上的部分也间隔插入导频，避免后面内插NaN
msg_pilot_aft_casd01_2(1, num_carrier * group_prl) = 1;   %补1是为了把不上的部分也间隔插入导频，避免后面内插NaN
msg_pilot_aft_casd10_1(1, num_carrier * group_prl) = 1;   %补1是为了把不上的部分也间隔插入导频，避免后面内插NaN
msg_pilot_aft_casd10_2(1, num_carrier * group_prl) = 1;   %补0是为了把不上的部分也间隔插入导频，避免后面内插NaN

msg_prl01_1 = transpose(reshape(msg_pilot_aft_casd01_1, num_carrier, group_prl));     %parallel signals
msg_prl01_2 = transpose(reshape(msg_pilot_aft_casd01_2, num_carrier, group_prl));     %parallel signals
msg_prl10_1 = transpose(reshape(msg_pilot_aft_casd10_1, num_carrier, group_prl));     %parallel signals
msg_prl10_2 = transpose(reshape(msg_pilot_aft_casd10_2, num_carrier, group_prl));     %parallel signals

PLoc_1 = 1: (pilotfrequency+1) :length(msg_pilot_aft_casd01_1);   %导频位置
DLoc_1 = setxor(1:length(msg_pilot_aft_casd01_1), PLoc_1);      %数据位置

%% OFDM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = sqrt(num_carrier);  %normalization factor
msg_idft10_1 = zeros(group_prl, num_carrier);
msg_idft10_2 = zeros(group_prl, num_carrier);
msg_idft01_1 = zeros(group_prl, num_carrier);
msg_idft01_2 = zeros(group_prl, num_carrier);

msg_ofdm10_1 = zeros(group_prl, num_carrier+len_cp);
msg_ofdm10_2 = zeros(group_prl, num_carrier+len_cp);
msg_ofdm01_1 = zeros(group_prl, num_carrier+len_cp);
msg_ofdm01_2 = zeros(group_prl, num_carrier+len_cp);

msg_ray0 = zeros(num_carrier+len_cp, group_prl);
msg_ray = zeros(group_prl, num_carrier+len_cp);

msg_nocp01_1 = zeros(group_prl, num_carrier);
msg_nocp01_2 = zeros(group_prl, num_carrier);
msg_nocp10_1 = zeros(group_prl, num_carrier);
msg_nocp10_2 = zeros(group_prl, num_carrier);

msg_nocp_ray = zeros(group_prl, num_carrier);
msg_nocp_est = zeros(group_prl, num_carrier);

msg_fft01_1 = zeros(group_prl, num_carrier);
msg_fft01_2 = zeros(group_prl, num_carrier);
msg_fft10_1 = zeros(group_prl, num_carrier);
msg_fft10_2 = zeros(group_prl, num_carrier);

MSG_PRL = zeros(group_prl, num_carrier);
MSG_NOCP_EST = zeros(group_prl, num_carrier);
H_RAY = zeros(group_prl, num_carrier);
% h_ray = zeros(group_prl, num_carrier+len_cp);
cp01_1 = zeros(group_prl, len_cp);  %initial
cp01_2 = zeros(group_prl, len_cp);  %initial
cp10_1 = zeros(group_prl, len_cp);  %initial
cp10_2 = zeros(group_prl, len_cp);  %initial

for j = 1 : group_prl
   msg_idft01_1(j,:) = ifft(msg_prl01_1(j,:))/N;    %idft
   msg_idft01_2(j,:) = ifft(msg_prl01_2(j,:))/N;    %idft
   msg_idft10_1(j,:) = ifft(msg_prl10_1(j,:))/N;    %idft   
   msg_idft10_2(j,:) = ifft(msg_prl10_2(j,:))/N;    %idft

   cp01_1(j,:) = msg_idft01_1(j,num_carrier-len_cp+1:num_carrier);   %cp, 12bit
   cp01_2(j,:) = msg_idft01_2(j,num_carrier-len_cp+1:num_carrier);   %cp, 12bit   
   cp10_1(j,:) = msg_idft10_1(j,num_carrier-len_cp+1:num_carrier);   %cp, 12bit
   cp10_2(j,:) = msg_idft10_2(j,num_carrier-len_cp+1:num_carrier);   %cp, 12bit
    
   msg_ofdm01_1(j,:) = cat(2,cp01_1(j,:),msg_idft01_1(j,:));    %add cp,concatenate matrix
   msg_ofdm01_2(j,:) = cat(2,cp01_2(j,:),msg_idft01_2(j,:));    %add cp,concatenate matrix
   msg_ofdm10_1(j,:) = cat(2,cp10_1(j,:),msg_idft10_1(j,:));    %add cp,concatenate matrix
   msg_ofdm10_2(j,:) = cat(2,cp10_2(j,:),msg_idft10_2(j,:));    %add cp,concatenate matrix
   
end
   %% MIMO_chan%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           
%     MIMO_chan = comm.MIMOChannel(...
%         'SampleRate', chanSRate,...
%         'MaximumDopplerShift', Doppler, ...
%         'PathDelays', PathDelays,...
%         'AveragePathGains', PathGains,...
%         'RandomStream', 'mt19937ar with seed',...
%         'Seed', 100,...
%         'SpatialCorrelation', false,...
%         'NumTransmitAntennas', numTx,...
%         'TransmitCorrelationMatrix', eye(numTx),...
%         'NumReceiveAntennas', numRx,...
%         'ReceiveCorrelationMatrix', eye(numRx),...
%         'PathGainsOutputPort', true,...
%         'NormalizePathGains', true,...
%         'NormalizeChannelOutputs', true);

msg_ofdm_re10_1 = reshape(transpose(msg_ofdm10_1),size(msg_ofdm10_1,1)*size(msg_ofdm10_1,2),1);   %每一列是一个天线的输入
msg_ofdm_re10_2 = reshape(transpose(msg_ofdm10_2),size(msg_ofdm10_2,1)*size(msg_ofdm10_2,2),1);   %每一列是一个天线的输入
msg_ofdm_re01_1 = reshape(transpose(msg_ofdm01_1),size(msg_ofdm01_1,1)*size(msg_ofdm01_1,2),1);   %每一列是一个天线的输入
msg_ofdm_re01_2 = reshape(transpose(msg_ofdm01_2),size(msg_ofdm01_2,1)*size(msg_ofdm01_2,2),1);   %每一列是一个天线的输入

%%%%%%%%%%%%%%%%% 4 Rayleigh represent MIMO%%%%%%%%%%%%%%%%%%    
rng(1213);    
Ray_chan11 = rayleighchan(Ts, Doppler, PathDelays, PathGains);    %rayleigh channel 
%rayleigh channel H11
rng(1263);    
Ray_chan12 = rayleighchan(Ts, Doppler, PathDelays, PathGains);     %rayleigh channel H12

rng(1923);    
Ray_chan21 = rayleighchan(Ts, Doppler, PathDelays, PathGains);     %rayleigh channel H21    

rng(109);    
Ray_chan22 = rayleighchan(Ts, Doppler, PathDelays, PathGains);     %rayleigh channel H22 

msg_ray11_01_1 = filter(Ray_chan11, msg_ofdm_re01_1);     %h11*x1,x1=0
msg_ray11_10_1 = filter(Ray_chan11, msg_ofdm_re10_1);     %h11*x1,x1=1

msg_ray12_01_2 = filter(Ray_chan12, msg_ofdm_re01_2);     %h12*x2,x2=1
msg_ray12_10_2 = filter(Ray_chan12, msg_ofdm_re10_2);     %h12*x2,x2=0

msg_MIMO01_1 = msg_ray11_01_1 + msg_ray12_01_2;      %y1 = h11*x1 + h12*x2,x1=0,x2=1
msg_MIMO10_1 = msg_ray11_10_1 + msg_ray12_10_2;      %y1 = h11*x1 + h12*x2,x1=1,x2=0

msg_ray21_01 = filter(Ray_chan21, msg_ofdm_re01_1);     %h21*x1,x1=0
msg_ray21_10 = filter(Ray_chan21, msg_ofdm_re10_1);     %h21*x1,x1=1

msg_ray22_01 = filter(Ray_chan22, msg_ofdm_re01_2);     %h22*x2,x2=1
msg_ray22_10 = filter(Ray_chan22, msg_ofdm_re10_2);     %h22*x2,x2=0

msg_MIMO01_2 = msg_ray21_01 + msg_ray22_01;      %y2 = h21*x1 + h22*x2,x1=0,x2=1
msg_MIMO10_2 = msg_ray21_10 + msg_ray22_10;      %y2 = h21*x1 + h22*x2,x1=1,x2=0

%     [msg_MIMO, MIMO_PATHGAINS] = step(MIMO_chan, msg_ofdm_re);   %pass MIMO, y = Hx
msg_awgn01_1 = awgn(msg_MIMO01_1,SNR_awgn,'measured',19980611);      %pass awgn,接收信号y = Hx + n
msg_awgn01_2 = awgn(msg_MIMO01_2,SNR_awgn,'measured',19980611);      %pass awgn,接收信号y = Hx + n
msg_awgn10_1 = awgn(msg_MIMO10_1,SNR_awgn,'measured',19980611);      %pass awgn,接收信号y = Hx + n
msg_awgn10_2 = awgn(msg_MIMO10_2,SNR_awgn,'measured',19980611);      %pass awgn,接收信号y = Hx + n

msg_MIMO_re01_1 = transpose(reshape(msg_awgn01_1, size(msg_awgn01_1,1)*size(msg_awgn01_1,2)/group_prl,group_prl));   %返回到并行子载波形式
msg_MIMO_re01_2 = transpose(reshape(msg_awgn01_2, size(msg_awgn01_2,1)*size(msg_awgn01_2,2)/group_prl,group_prl));   %返回到并行子载波形式
msg_MIMO_re10_1 = transpose(reshape(msg_awgn10_1, size(msg_awgn10_1,1)*size(msg_awgn10_1,2)/group_prl,group_prl));   %返回到并行子载波形式
msg_MIMO_re10_2 = transpose(reshape(msg_awgn10_2, size(msg_awgn10_2,1)*size(msg_awgn10_2,2)/group_prl,group_prl));   %返回到并行子载波形式

for j = 1 : group_prl       
   %% DFT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   msg_nocp01_1(j,:) = msg_MIMO_re01_1(j,len_cp+1:num_carrier+len_cp);    %wipe off cp
   msg_nocp01_2(j,:) = msg_MIMO_re01_2(j,len_cp+1:num_carrier+len_cp);    %wipe off cp   
   msg_nocp10_1(j,:) = msg_MIMO_re10_1(j,len_cp+1:num_carrier+len_cp);    %wipe off cp
   msg_nocp10_2(j,:) = msg_MIMO_re10_2(j,len_cp+1:num_carrier+len_cp);    %wipe off cp

   msg_fft01_1(j,:) = N*fft(msg_nocp01_1(j,:));     %idft
   msg_fft01_2(j,:) = N*fft(msg_nocp01_2(j,:));     %idft
   msg_fft10_1(j,:) = N*fft(msg_nocp10_1(j,:));     %idft
   msg_fft10_2(j,:) = N*fft(msg_nocp10_2(j,:));     %idft


end

%% 并串变换%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

msg_casd01_1 = transpose(reshape(transpose(msg_fft01_1),size(msg_fft01_1,1)*size(msg_fft01_1,2),1));  %y1, x1 = 0,x2 = 1接收，h12
msg_casd01_2 = transpose(reshape(transpose(msg_fft01_2),size(msg_fft01_2,1)*size(msg_fft01_2,2),1));  %y2, x1 = 0,x2 = 1接收, h22
msg_casd10_1 = transpose(reshape(transpose(msg_fft10_1),size(msg_fft10_1,1)*size(msg_fft10_1,2),1));  %y1, x1 = 1,x2 = 0接收，h11
msg_casd10_2 = transpose(reshape(transpose(msg_fft10_2),size(msg_fft10_2,1)*size(msg_fft10_2,2),1));  %y2, x1 = 1,x2 = 0接收，h21

%% 信道估计%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%LS channel estimate%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h01_1 = msg_casd01_1(PLoc_1);     %notice that pilot is 1
h01_2 = msg_casd01_2(PLoc_1);     %notice that pilot is 1
h10_1 = msg_casd10_1(PLoc_1);     %notice that pilot is 1
h10_2 = msg_casd10_2(PLoc_1);     %notice that pilot is 1

msg_casd_len = length(msg_casd01_1);

h_est01_1 = interp1(PLoc_1, h01_1, 1:msg_casd_len);
h_est01_2 = interp1(PLoc_1, h01_2, 1:msg_casd_len);
h_est10_1 = interp1(PLoc_1, h10_1, 1:msg_casd_len);
h_est10_2 = interp1(PLoc_1, h10_2, 1:msg_casd_len);

h_est01_1(msg_casd_len - rest_num + 1 :msg_casd_len) = [];   %去掉补的0
h_est01_2(msg_casd_len - rest_num + 1 :msg_casd_len) = [];   %去掉补的0
h_est10_1(msg_casd_len - rest_num + 1 :msg_casd_len) = [];   %去掉补的0
h_est10_2(msg_casd_len - rest_num + 1 :msg_casd_len) = [];   %去掉补的0

PLoc = 1: (pilotfrequency+1) :length(h_est01_1); 
DLoc = setxor(1:length(h_est01_1), PLoc);      %数据位置
h_ray01_1 = h_est01_1(DLoc);
h_ray01_2 = h_est01_2(DLoc);
h_ray10_1 = h_est10_1(DLoc);
h_ray10_2 = h_est10_2(DLoc);

h11 = h_ray10_1;
h12 = h_ray01_1;
h21 = h_ray10_2;
h22 = h_ray01_2;    %信道估计值

msg_casd01_1(msg_casd_len - rest_num + 1 :msg_casd_len) = [];     %去掉串并变换补的0
msg_casd01_1 = msg_casd01_1(DLoc);
msg_casd01_2(msg_casd_len - rest_num + 1 :msg_casd_len) = [];     %去掉串并变换补的0
msg_casd01_2 = msg_casd01_2(DLoc);
msg_casd10_1(msg_casd_len - rest_num + 1 :msg_casd_len) = [];     %去掉串并变换补的0
msg_casd10_1 = msg_casd10_1(DLoc);
msg_casd10_2(msg_casd_len - rest_num + 1 :msg_casd_len) = [];     %去掉串并变换补的0
msg_casd10_2 = msg_casd10_2(DLoc);

msg_casd10 = cat(2, msg_casd10_1, msg_casd10_2);
msg_casd10_re = reshape(msg_casd10, length(msg_casd10)/2,2);  %转化问天线输出形式
%%  MIMO接收%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = [h11,h12;h21,h22];      %信道矩阵
x_base = (1/sqrt(2))*[1+1i, 1+1i];
x = zeros(16, 2);   %qpsk所有可能，16种，用x_base生成
for j = 1: 4
    for i = 0:3
        x(4*i+j, 1) = x_base(1) * exp(i*1i*pi/2);
        x(4*i+j, 2) = x_base(2) * exp(j*1i*pi/2);
    end        
end

y = zeros(16, 2);   %y is Hx, 16种
flag = zeros(1,16);
msg_MIMO_rx = zeros(size(msg_casd10_re,1)*size(msg_casd10_re,2)/2, 2);    %ML后的接收信号
for i = 1 : size(msg_casd10_re,1)*size(msg_casd10_re,2)/2
    for j = 1:16
        y(j,1) = h11(i)*x(j,1) + h12(i)*x(j,2);
        y(j,2) = h21(i)*x(j,1) + h22(i)*x(j,2);
        flag(j) = norm(msg_casd10_re(i,:)-y(j,:));   %in order to compare 求范数=-=立个flag,y is Hx
    end
    [min_num, min_pos] = min(flag);     %找到最小flag的索引
    msg_MIMO_rx(i,:) = x(min_pos,:);    %ML complete!!!
end
msg_MIMO_rx_re = reshape(msg_MIMO_rx, 1, size(msg_MIMO_rx,1)*size(msg_MIMO_rx,2));  %转化为行向量

    

%% QPSK解调%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg_gry_demod = pskdemod(msg_MIMO_rx_re,M, pi/4);      %qpsk demodulate
[dummy, graydecode] = sort(grayencode); 
graydecode = graydecode - 1;
msg_demod = graydecode(msg_gry_demod+1)';
msg_demod_bi = de2bi(msg_demod,k)'; 
msg_demod_bi = msg_demod_bi(:);     %gray->normal

%% 解行列交织%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% msg_deintrlv_re = transpose(reshape(msg_demod_bi,row_intrlv*column_intrlv,group_num));      %row_intrlv*column_intrlvbits as a group again
% msg_deintrlv = zeros(group_num, row_intrlv*column_intrlv);
% for k = 1:group_num
%     msg_deintrlv(k,:) = matdeintrlv(msg_deintrlv_re(k,:),row_intrlv, column_intrlv);   %relieve interleave
% end 
% msg_deintrlv_code = reshape(transpose(msg_deintrlv),1,size(msg_deintrlv,1)*size(msg_deintrlv,2));   %code after deinterleaving    
% deintrlv_len = length(msg_deintrlv_code);
% msg_deintrlv_code(deintrlv_len - rest_num1 + 1 : deintrlv_len) = [];     %eliminate 0


%% 解随机交织%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg_deintrlv_re = transpose(reshape(msg_demod_bi,intrl_frm, intrl_rand_ts));
msg_deintrlv = zeros(intrl_rand_ts, intrl_frm);
for k = 1:intrl_rand_ts
    msg_deintrlv(k, :) = randdeintrlv(msg_deintrlv_re(k, :), k);
end
msg_deintrlv_code = reshape(transpose(msg_deintrlv),1,size(msg_deintrlv,1)*size(msg_deintrlv,2));   %code after deinterleaving    
deintrl_len = length(msg_deintrlv_code);
msg_deintrlv_code(deintrl_len - rest_num1 + 1 : deintrl_len) = [];     %去掉0

% msg_deintrlv_code = msg_demod_bi'; %test

%% 解卷积码%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg_rx = vitdec(msg_deintrlv_code, trellis, tblen, 'cont', 'hard');     %'hard' judgment
msg_rx(1:35) = [];  %delete traceback daley



%% 解Turbo码%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turbo_decoder = comm.TurboDecoder(...
%     'TrellisStructure', poly2trellis(4, [13 15 17], 13),...
%     'InterleaverIndices', intrlvrIndices,...
%     'NumIterations', 6);
% msg_rx = step(Turbo_decoder, transpose(msg_deintrlv_code));     %Turbo decode
% msg_rx = transpose(msg_rx);

%% 计算BER&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
msg(num_bit-34: num_bit) = [];      %for convolutional code
% msg_rx = msg_deintrlv_code;
msg_wrong = abs(msg_rx - msg);      %for turbo code
BER = sum(msg_wrong)/num_bit;