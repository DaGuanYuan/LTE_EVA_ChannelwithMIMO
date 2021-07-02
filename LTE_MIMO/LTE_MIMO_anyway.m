%LTE_7:交织器和译码器的选择研究
%% 初始条件和信号生成%%%%%%%%%%%%%%%%%%%%%%%%%%
num_bit = 1.28e5;
Nsamp = 4;
M = 4;
k = log2(M);
Rb = 20*(10^6);    %bit rate 
Rs = Rb / k;
Tb = 1/Rb;   
Ts = 1/Rs;
rng('default');
rng(67);    %seed = 67
msg = randi([0,1],1,num_bit);  %generate messege
SNR_awgn  = 1000000;%awgn snr
%%MIMO%%
numTx = 2;
numRx = 2;
chanSRate = Rs;
% PathDelays = [0,30,150,310,370,710,1090,1730,2510]*1e-9 ;
% PathGains =  [0,-1.5,-1.4,-3.6,-0.6,-9.1,-7.0,-12.0,-16.9];
PathDelays = [0];
PathGains = [0];
velocity = 60/3.6;      %velocity
f = 2*10^9;     %载频
c = 3*10^8;     %light speed
Doppler = velocity*f/c;     %Maxinum Doppler shift




%% %%%%%%%卷积编码%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constlen = 7;   %constraint length
% codegen = [171, 173];   % polynomials generator
% tblen = 5*constlen;     %traceback length
% trellis = poly2trellis(constlen, codegen);
% dspec = distspec(trellis, 7);
% msg_conv = convenc(msg, trellis);   %convolutional encode
% msg_turbo_cascade = msg_conv;%Here msg_conv_re could be both turbo encode and convolution encode
% G = length(msg_conv)/length(msg);    %code generator gain of the convolutional encode


%% Turbo编码 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frmlen = num_bit;     %The frame length ->may not be useful
intrlvrIndices = randperm(num_bit);
Turbo_encoder = comm.TurboEncoder(...
    'TrellisStructure', poly2trellis(4, [13 15 17], 13),...
    'InterleaverIndices', intrlvrIndices);  %This is the Turbo encoder
msg_turbo = step(Turbo_encoder, transpose(msg));
msg_turbo_cascade = transpose(msg_turbo);
msg_conv = msg_turbo_cascade;  %Here msg_conv_re could be both turbo encode and convolution encode
G = length(msg_turbo_cascade)/length(msg);  %code generator gain of the turbo encode


%% 行列交织 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row_intrlv = 1;     %interleave depth
column_intrlv = 128;   %interleave width
group_num = fix(length(msg_conv)/(row_intrlv*column_intrlv));  %row_intrlv*column_intrlv bits as a group;group number
more_num1 = length(msg_conv) - row_intrlv*column_intrlv*group_num;
rest_num1 = row_intrlv*column_intrlv* (group_num + 1) - length(msg_conv);
group_num = group_num + 1;
msg_conv(1, group_num*row_intrlv*column_intrlv) = 0;    %补0
msg_conv_re = transpose(reshape(msg_conv,row_intrlv*column_intrlv,group_num));     %group by group
msg_intrlv = zeros(group_num,row_intrlv*column_intrlv);     %initial->distribute caches
for i = 1:group_num
    msg_intrlv(i,:) = matintrlv(msg_conv_re(i,:), row_intrlv, column_intrlv);    %100*50 matrix interleave
end
msg_intrlv_cascade = transpose(reshape(transpose(msg_intrlv), 1, size(msg_intrlv,1)*size(msg_intrlv,2)));   %interleave completed!!!!!!!!!!    



%% 随机交织%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intrl_frm = 1000;
% intrl_rand_ts =fix(length(msg_conv) / intrl_frm);   %全都补零部分情况了！！！！
% more_num1 = length(msg_conv) - intrl_frm*intrl_rand_ts;
% rest_num1 = intrl_frm*(intrl_rand_ts + 1) - length(msg_conv);  %补的0的个数
% intrl_rand_ts = intrl_rand_ts + 1;
% msg_turbo_cascade(1, intrl_frm*(intrl_rand_ts)) = 0;  %无论如何都补0，反正都去掉
% msg_conv_re = transpose(reshape(msg_turbo_cascade, intrl_frm, intrl_rand_ts));%here msg_conv_re could be both turbo encode and convolution encode
% msg_intrlv = zeros(intrl_rand_ts , intrl_frm);
% for i = 1:intrl_rand_ts
%     msg_intrlv(i,:) = randintrlv(msg_conv_re(i,:), i) ;     %i is the random seed, 我真是个天才=-=这样每一个的随机过程都不一样
% end
% msg_intrlv_cascade = transpose(reshape(transpose(msg_intrlv), 1, size(msg_intrlv,1)*size(msg_intrlv,2)));  %interleave completed!!!!!!!!!!!
% 


%% QPSK调制%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg_syb = transpose(reshape(msg_intrlv_cascade, size(msg_intrlv_cascade, 2)*k, size(msg_intrlv_cascade, 1)/k));     %IQ, get symbols
msg_de = bi2de(msg_syb);    %decimal symbols
grayencode = bitxor(0:M-1, floor((0:M-1)/2));   %generate gray list by order,[0,1,3,2]
msg_gry_enc = grayencode(msg_de+1);     %gray encode:0->0(The first),1->1(The second)，2->3(The third),3->2(The forth)，so by order,msg_enc+1
msg_qpsk = pskmod(msg_gry_enc,M, pi/4);    %QPSK modulation

%% 串并变换%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_carrier = 128;     %子载波数量
len_cp = 12;    %length of cp


%% 串并变换%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
group_prl = fix(length(msg_qpsk)/num_carrier);   %number of the groups which are divided in serial-to-parallel conversion
more_num3 = length(msg_qpsk) - num_carrier*group_prl; %the same as interleaving
rest_num3 = num_carrier*(group_prl + 1) - length(msg_qpsk);  %第二次补0个数
rest_num = rest_num2 + rest_num3 ;   %总补零的个数
group_prl = group_prl + 1;
msg_qpsk(1, num_carrier * group_prl) = 1;   %补1是为了把不上的部分也间隔插入导频，避免后面内插NaN

msg_prl = transpose(reshape(msg_qpsk, num_carrier, group_prl));     %parallel signals
PLoc_1 = 1: (pilotfrequency+1) :length(msg_qpsk);   %导频位置
DLoc_1 = setxor(1:length(msg_qpsk), PLoc_1);      %数据位置

%% OFDM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = sqrt(num_carrier);  %normalization factor
msg_idft = zeros(group_prl, num_carrier);
msg_ofdm = zeros(group_prl, num_carrier+len_cp);
msg_ray0 = zeros(num_carrier+len_cp, group_prl);
msg_ray = zeros(group_prl, num_carrier+len_cp);
msg_nocp = zeros(group_prl, num_carrier);
msg_nocp_ray = zeros(group_prl, num_carrier);
msg_nocp_est = zeros(group_prl, num_carrier);
msg_fft = zeros(group_prl, num_carrier);
MSG_PRL = zeros(group_prl, num_carrier);
MSG_NOCP_EST = zeros(group_prl, num_carrier);
H_RAY = zeros(group_prl, num_carrier);
% h_ray = zeros(group_prl, num_carrier+len_cp);
cp = zeros(group_prl, 12);  %initial
for j = 1 : group_prl
   msg_idft(j,:) = ifft(msg_prl(j,:))/N;    %idft
   cp(j,:) = msg_idft(j,117:num_carrier);   %cp, 12bit
   msg_ofdm(j,:) = cat(2,cp(j,:),msg_idft(j,:));    %add cp,concatenate matrix
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

msg_ofdm_re = reshape(transpose(msg_ofdm),size(msg_ofdm,1)*size(msg_ofdm,2)/numTx,numTx);   %每一列是一个天线的输入
%%%%%%%%%%%%%%%%% 4 Rayleigh represent MIMO%%%%%%%%%%%%%%%%%%    
Ray_chan11 = comm.RayleighChannel(...
    'SampleRate', chanSRate,...
    'PathDelays', PathDelays,...
    'AveragePathGains', PathGains,...
    'NormalizePathGains', true,...
    'MaximumDopplerShift', Doppler,...
    'RandomStream','mt19937ar with seed', ...
    'Seed', 50,...
    'PathGainsOutputPort', true);     %rayleigh channel H11

Ray_chan12 = comm.RayleighChannel(...
    'SampleRate', chanSRate,...
    'PathDelays', PathDelays,...
    'AveragePathGains', PathGains,...
    'NormalizePathGains', true,...
    'MaximumDopplerShift', Doppler,...
    'RandomStream','mt19937ar with seed', ...
    'Seed', 100,...
    'PathGainsOutputPort', true);     %rayleigh channel H12

Ray_chan21 = comm.RayleighChannel(...
    'SampleRate', chanSRate,...
    'PathDelays', PathDelays,...
    'AveragePathGains', PathGains,...
    'NormalizePathGains', true,...
    'MaximumDopplerShift', Doppler,...
    'RandomStream','mt19937ar with seed', ...    
    'Seed', 150,...
    'PathGainsOutputPort', true);     %rayleigh channel H21    

Ray_chan22 = comm.RayleighChannel(...
    'SampleRate', chanSRate,...
    'PathDelays', PathDelays,...
    'AveragePathGains', PathGains,...
    'NormalizePathGains', true,...
    'MaximumDopplerShift', Doppler,...
    'RandomStream','mt19937ar with seed', ...
    'Seed', 200,...
    'PathGainsOutputPort', true);     %rayleigh channel H22 

msg_MIMO = zeros(size(msg_ofdm_re,1) ,2);
[msg_ray11,h11] = step(Ray_chan11, msg_ofdm_re(:,1));     %h11*x1
[msg_ray12,h12] = step(Ray_chan12, msg_ofdm_re(:,2));     %h12*x2
msg_MIMO(:,1) = msg_ray11 + msg_ray12;      %y1 = h11*x1 + h12*x2
[msg_ray21,h21] = step(Ray_chan21, msg_ofdm_re(:, 1));     %h21*x1
[msg_ray22,h22] = step(Ray_chan22, msg_ofdm_re(:, 2));     %h22*x2
msg_MIMO(:,2) = msg_ray21 + msg_ray22;      %y2 = h21*x1 + h22*x2
h1 = cat(1,h11,h12);
h2 = cat(1,h21,h22);
h1 = transpose(reshape(h1, length(h1)/group_prl,group_prl));
h2 = transpose(reshape(h2, length(h2)/group_prl,group_prl));

% h11 = transpose(reshape(h11, length(h11)/group_prl,group_prl));
% h12 = transpose(reshape(h12, length(h12)/group_prl,group_prl));
% h21 = transpose(reshape(h21, length(h21)/group_prl,group_prl));
% h22 = transpose(reshape(h22, length(h22)/group_prl,group_prl));


%     [msg_MIMO, MIMO_PATHGAINS] = step(MIMO_chan, msg_ofdm_re);   %pass MIMO, y = Hx
msg_awgn = awgn(msg_MIMO,SNR_awgn,'measured',19980611);      %pass awgn,接收信号y = Hx + n
msg_MIMO_re = transpose(reshape(msg_awgn, size(msg_awgn,1)*size(msg_awgn,2)/group_prl,group_prl));   %返回到并行子载波形式

for j = 1 : group_prl       
   %% DFT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   msg_nocp(j,:) = msg_MIMO_re(j,13:num_carrier+len_cp);    %wipe off cp
   h1_re(j,:) = h1(j,13:num_carrier+len_cp);
   h2_re(j,:) = h2(j,13:num_carrier+len_cp);
   
   
   msg_fft(j,:) = N*fft(msg_nocp(j,:));     %idft

end

%% 并串变换%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg_casd = reshape(transpose(msg_fft),size(msg_fft,1)*size(msg_fft,2)/2,2);   %按照接收mimo信号
h1_re = reshape(transpose(h1_re),size(h1_re,1)*size(h1_re,2)/2,2);
h2_re = reshape(transpose(h2_re),size(h2_re,1)*size(h2_re,2)/2,2);
h11_re = transpose(h1_re(:,1));
h12_re = transpose(h1_re(:,2));
h21_re = transpose(h2_re(:,1));
h22_re = transpose(h2_re(:,2));

%% 信道估计%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%LS channel estimate%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
msg_MIMO_rx = zeros(size(msg_casd,1)*size(msg_casd,2)/2, 2);    %ML后的接收信号
for i = 1 : size(msg_casd,1)*size(msg_casd,2)/2
    for j = 1:16
        y(j,1) = h11_re(i)*x(j,1) + h12_re(i)*x(j,2);
        y(j,2) = h21_re(i)*x(j,1) + h22_re(i)*x(j,2);
        flag(j) = norm(msg_casd(i,:)-y(j,:));   %in order to compare 求范数=-=立个flag,y is Hx
    end
    [min_num, min_pos] = min(flag);     %找到最小flag的索引
    msg_MIMO_rx(i,:) = x(min_pos,:);    %ML complete!!!
end
msg_MIMO_rx_re = reshape(msg_MIMO_rx, 1, size(msg_MIMO_rx,1)*size(msg_MIMO_rx,2));  %转化为行向量
msg_MIMO_len = length(msg_MIMO_rx_re);
msg_MIMO_rx_re(msg_MIMO_len - rest_num3 +1 : msg_MIMO_len) = [];     %去掉补的0
    

%% QPSK解调%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg_gry_demod = pskdemod(msg_MIMO_rx_re,M, pi/4);      %qpsk demodulate
[dummy, graydecode] = sort(grayencode); 
graydecode = graydecode - 1;
msg_demod = graydecode(msg_gry_demod+1)';
msg_demod_bi = de2bi(msg_demod,k)'; 
msg_demod_bi = msg_demod_bi(:);     %gray->normal

%% 解行列交织%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg_deintrlv_re = transpose(reshape(msg_demod_bi,row_intrlv*column_intrlv,group_num));      %row_intrlv*column_intrlvbits as a group again
msg_deintrlv = zeros(group_num, row_intrlv*column_intrlv);
for k = 1:group_num
    msg_deintrlv(k,:) = matdeintrlv(msg_deintrlv_re(k,:),row_intrlv, column_intrlv);   %relieve interleave
end 
msg_deintrlv_code = reshape(transpose(msg_deintrlv),1,size(msg_deintrlv,1)*size(msg_deintrlv,2));   %code after deinterleaving    
deintrlv_len = length(msg_deintrlv_code);
msg_deintrlv_code(deintrlv_len - rest_num1 + 1 : deintrlv_len) = [];     %eliminate 0


%% 解随机交织%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% msg_deintrlv_re = transpose(reshape(msg_demod_bi,intrl_frm, intrl_rand_ts));
% msg_deintrlv = zeros(intrl_rand_ts, intrl_frm);
% for k = 1:intrl_rand_ts
%     msg_deintrlv(k, :) = randdeintrlv(msg_deintrlv_re(k, :), k);
% end
% msg_deintrlv_code = reshape(transpose(msg_deintrlv),1,size(msg_deintrlv,1)*size(msg_deintrlv,2));   %code after deinterleaving    
% deintrl_len = length(msg_deintrlv_code);
% msg_deintrlv_code(deintrl_len - rest_num1 + 1 : deintrl_len) = [];     %去掉0



%% 解卷积码%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% msg_rx = vitdec(msg_deintrlv_code, trellis, tblen, 'cont', 'hard');     %'hard' judgment
% msg_rx(1:35) = [];  %delete traceback daley



%% 解Turbo码%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Turbo_decoder = comm.TurboDecoder(...
    'TrellisStructure', poly2trellis(4, [13 15 17], 13),...
    'InterleaverIndices', intrlvrIndices,...
    'NumIterations', 6);
msg_rx = step(Turbo_decoder, transpose(msg_deintrlv_code));     %Turbo decode
msg_rx = transpose(msg_rx);

%% 计算BER&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% msg(num_bit-34: num_bit) = [];      %for convolutional code
msg_wrong = abs(msg_rx - msg);      %for turbo code
BER = sum(msg_wrong)/num_bit;


