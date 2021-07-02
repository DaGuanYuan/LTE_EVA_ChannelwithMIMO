%LTE_7:交织器和译码器的选择研究
%%%%%%%%%%%初始条件和信号生成%%%%%%%%%%%%%%%%%%%%%%%%%%
num_bit = 1.28e5;
Nsamp = 4;
M = 4;
k = log2(M);
Rb = 200000;    %bit rate 
Rs = Rb / k;
Tb = 1/Rb;   
Ts = 1/Rs;
rng('default');
rng(67);    %seed = 67
msg = randi([0,1],1,num_bit);  %generate messege
SNR_awgn  = 0;%awgn snr
snr = SNR_awgn;

% %%%%%%%%%卷积编码%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% constlen = 7;   %constraint length
% codegen = [171, 173];   % polynomials generator
% tblen = 5*constlen;     %traceback length
% trellis = poly2trellis(constlen, codegen);
% dspec = distspec(trellis, 7);
% msg_conv = convenc(msg, trellis);   %convolutional encode
% msg_turbo_cascade = msg_conv;%Here msg_conv_re could be both turbo encode and convolution encode
% G = length(msg_conv)/length(msg);    %code generator gain of the convolutional encode


%%%%%%%%%%%Turbo编码%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frmlen = num_bit;     %The frame length ->may not be useful
intrlvrIndices = randperm(num_bit);
Turbo_encoder = comm.TurboEncoder(...
    'TrellisStructure', poly2trellis(4, [13 15 17], 13),...
    'InterleaverIndices', intrlvrIndices);  %This is the Turbo encoder
msg_turbo = step(Turbo_encoder, transpose(msg));
msg_turbo_cascade = transpose(msg_turbo);
msg_conv = msg_turbo_cascade;  %Here msg_conv_re could be both turbo encode and convolution encode
G = length(msg_turbo_cascade)/length(msg);  %code generator gain of the turbo encode


%%%%%%%%%%%%行列交织%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
row_intrlv = 1;     %interleave depth
column_intrlv = 128;   %interleave width
group_num = fix(length(msg_conv)/(row_intrlv*column_intrlv));  %row_intrlv*column_intrlv bits as a group;group number
more_num1 = length(msg_conv) - row_intrlv*column_intrlv*group_num;
rest_num1 = row_intrlv*column_intrlv* (group_num + 1) - length(msg_conv);
if more_num1 ~= 0
    group_num = group_num + 1;
    msg_conv(1, group_num*row_intrlv*column_intrlv) = 0;    %补0
    msg_conv_re = transpose(reshape(msg_conv,row_intrlv*column_intrlv,group_num));     %group by group
    msg_intrlv = zeros(group_num,row_intrlv*column_intrlv);     %initial->distribute caches
    for i = 1:group_num
        msg_intrlv(i,:) = matintrlv(msg_conv_re(i,:), row_intrlv, column_intrlv);    %100*50 matrix interleave
    end
    msg_intrlv_cascade = transpose(reshape(transpose(msg_intrlv), 1, size(msg_intrlv,1)*size(msg_intrlv,2)));   %interleave completed!!!!!!!!!!    
else
    msg_conv_re = transpose(reshape(msg_conv,row_intrlv*column_intrlv,group_num));     %group by group
    msg_intrlv = zeros(group_num,row_intrlv*column_intrlv);     %initial->distribute caches
    for i = 1:group_num
        msg_intrlv(i,:) = matintrlv(msg_conv_re(i,:), row_intrlv, column_intrlv);    %100*50 matrix interleave
    end
    msg_intrlv_cascade = transpose(reshape(transpose(msg_intrlv), 1, size(msg_intrlv,1)*size(msg_intrlv,2)));   %interleave completed!!!!!!!!!!    
end


% %%%%%%%%%%%%随机交织%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% intrl_frm = 1000;
% intrl_rand_ts =fix(length(msg_conv) / intrl_frm);   %分两块，除的进的和剩下的分开交织
% more_num1 = length(msg_conv) - intrl_frm*intrl_rand_ts;
% rest_num1 = intrl_frm*(intrl_rand_ts + 1) - length(msg_conv);  %补的0的个数
% if more_num1~=0
%     intrl_rand_ts = intrl_rand_ts + 1;
%     msg_turbo_cascade(1, intrl_frm*(intrl_rand_ts)) = 0;  %补0
%     msg_conv_re = transpose(reshape(msg_turbo_cascade, intrl_frm, intrl_rand_ts));%here msg_conv_re could be both turbo encode and convolution encode
%     msg_intrlv = zeros(intrl_rand_ts , intrl_frm);
%     for i = 1:intrl_rand_ts
%         msg_intrlv(i,:) = randintrlv(msg_conv_re(i,:), i) ;     %i is the random seed, 我真是个天才=-=这样每一个的随机过程都不一样
%     end
%     msg_intrlv_cascade = transpose(reshape(transpose(msg_intrlv), 1, size(msg_intrlv,1)*size(msg_intrlv,2)));  %interleave completed!!!!!!!!!!!
% 
% else
%     msg_conv_re = transpose(reshape(msg_turbo_cascade, (length(msg_conv))/intrl_rand_ts, intrl_rand_ts));  %here msg_conv_re could be both turbo encode and convolution encode
%     msg_intrlv = zeros(intrl_rand_ts, (length(msg_conv))/intrl_rand_ts);
%     for i = 1:intrl_rand_ts
%         msg_intrlv(i,:) = randintrlv(msg_conv_re(i,:), i) ;     %i is the random seed, 我真是个天才=-=这样每一个的随机过程都不一样
%     end
%     msg_intrlv_cascade = transpose(reshape(transpose(msg_intrlv), 1, length(msg_conv)));  %interleave completed!!!!!!!!!!!
% end
% 

%%%%%%%%%%%QPSK调制%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg_syb = transpose(reshape(msg_intrlv_cascade, size(msg_intrlv_cascade, 2)*k, size(msg_intrlv_cascade, 1)/k));     %IQ, get symbols
msg_de = bi2de(msg_syb);    %decimal symbols
grayencode = bitxor(0:M-1, floor((0:M-1)/2));   %generate gray list by order,[0,1,3,2]
msg_gry_enc = grayencode(msg_de+1);     %gray encode:0->0(The first),1->1(The second)，2->3(The third),3->2(The forth)，so by order,msg_enc+1
msg_qpsk = pskmod(msg_gry_enc,M, pi/4);    %QPSK modulation

%%%%%%%%%%%%串并变换%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_carrier = 128;     %子载波数量
len_cp = 12;    %length of cp
group_prl = fix(length(msg_qpsk)/num_carrier);   %number of the groups which are divided in serial-to-parallel conversion
more_num2 = length(msg_qpsk) - num_carrier*group_prl; %the same as interleaving
rest_num2 = num_carrier*(group_prl + 1) - length(msg_qpsk);  %补0个数
if more_num2 ~= 0
    group_prl = group_prl + 1;
    msg_qpsk(1, num_carrier * group_prl) = 0;   %补0
    msg_prl = transpose(reshape(msg_qpsk, num_carrier, group_prl));     %parallel signals

    %%%%%%%%%%%%OFDM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = sqrt(num_carrier);  %normalization factor
    msg_idft = zeros(group_prl, num_carrier);
    msg_ofdm = zeros(group_prl, num_carrier+len_cp);
    msg_ray0 = zeros(num_carrier+len_cp, group_prl);
    msg_ray = zeros(group_prl, num_carrier+len_cp);
    msg_awgn = zeros(group_prl, num_carrier+len_cp);
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

       %%%%%%%%%%%%%%%Rayleigh_chan%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    rayleigh_chan = comm.RayleighChannel(...
    %        'SampleRate', Rs,...
    %        'PathDelays', [0, 1.5e-5],...
    %        'AveragePathGains',[0, -2],...
    %        'NormalizePathGains',true,...
    %        'MaximumDopplerShift',100,...
    %        'PathGainsOutputPort',true);%rayleigh channel 
    %    [msg_ray0(:,j),PathGains1] = rayleigh_chan(transpose(msg_ofdm(j,:)));     %signal after rayleigh channel   
    %    msg_ray(j,:) = transpose(msg_ray0(:,j));

       rayleigh_chan = rayleighchan(Ts, 100, [0, 1.5e-5], [0, -2]);     %rayleigh channel
       rayleigh_chan.StorePathGains = 1;
       msg_ray(j,:) = filter(rayleigh_chan, msg_ofdm(j, :));    %signal after rayleigh channel
       %%%%%%%%%%%%%%%AWGN_chan%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        msg_awgn(j,:) = awgn(msg_ofdm(j,:),SNR_awgn,'measured');     %not passing rayleigh
       msg_awgn(j,:) = awgn(msg_ray(j,:),SNR_awgn,'measured',19980611);      %after passing rayleigh
%        msg_awgn(j,:) = msg_ofdm(j,:);    %test, no noise

       %%%%%%%%%%%%%%%%%%test_no awgn%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    msg_nocp(j,:) = msg_ray(j,13:num_carrier+len_cp);    %wipe off cp
    %    msg_fft(j,:) = N*fft(msg_nocp(j,:));     %idft
       %%%%%%%%%%%%%%%DFT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       msg_nocp(j,:) = msg_awgn(j,13:num_carrier+len_cp);    %wipe off cp
       msg_fft(j,:) = N*fft(msg_nocp(j,:));     %idft

          %%%%%%%%%%%%%%%%%%理想信道估计%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       msg_nocp_ray(j,:) = msg_ray(j,13:num_carrier+len_cp);   %only rayleigh
       msg_nocp_est(j,:) = N*fft(msg_nocp_ray(j,:));

    end

    h_ray = msg_nocp_est./msg_prl;
    %%%%%%%%%%%%%%%%%%并串变换%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    msg_casd = reshape(transpose(msg_fft),1,size(msg_fft,1)*size(msg_fft,2));
    msg_casd_len = length(msg_casd);
    msg_casd(msg_casd_len - rest_num2 +1 : msg_casd_len) = [];     %去掉补的0
    
    % H_RAY_casd = reshape(transpose(H_RAY),1,num_bit);
    h_ray_casd = reshape(transpose(h_ray),1,size(h_ray,1)*size(h_ray,2));
    h_ray_casd_len = length(h_ray_casd);
    h_ray_casd(h_ray_casd_len - rest_num2 +1 : h_ray_casd_len) = [];     %去掉补的0

    %%%%%%%%%%%%%%%%%%%信道均衡(f)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    msg_casd_zf = msg_casd./h_ray_casd;
%     msg_casd_zf_len = length(msg_casd_zf);
%     msg_casd_zf(msg_casd_zf_len - rest_num2 +1 : msg_casd_zf_len) = [];
    
else
    msg_prl = transpose(reshape(msg_qpsk, num_carrier, group_prl));     %parallel signals

    %%%%%%%%%%%%OFDM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    N = sqrt(num_carrier);  %normalization factor
    msg_idft = zeros(group_prl, num_carrier);
    msg_ofdm = zeros(group_prl, num_carrier+len_cp);
    msg_ray0 = zeros(num_carrier+len_cp, group_prl);
    msg_ray = zeros(group_prl, num_carrier+len_cp);
    msg_awgn = zeros(group_prl, num_carrier+len_cp);
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

       %%%%%%%%%%%%%%%Rayleigh_chan%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    rayleigh_chan = comm.RayleighChannel(...
    %        'SampleRate', Rs,...
    %        'PathDelays', [0, 1.5e-5],...
    %        'AveragePathGains',[0, -2],...
    %        'NormalizePathGains',true,...
    %        'MaximumDopplerShift',100,...
    %        'PathGainsOutputPort',true);%rayleigh channel 
    %    [msg_ray0(:,j),PathGains1] = rayleigh_chan(transpose(msg_ofdm(j,:)));     %signal after rayleigh channel   
    %    msg_ray(j,:) = transpose(msg_ray0(:,j));

       rayleigh_chan = rayleighchan(Ts, 100, [0, 1.5e-5], [0, -2]);     %rayleigh channel
       rayleigh_chan.StorePathGains = 1;
       msg_ray(j,:) = filter(rayleigh_chan, msg_ofdm(j, :));    %signal after rayleigh channel
       %%%%%%%%%%%%%%%AWGN_chan%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    msg_awgn(j,:) = awgn(msg_ofdm(j,:),SNR_awgn,'measured');     %not passing rayleigh
       msg_awgn(j,:) = awgn(msg_ray(j,:),SNR_awgn,'measured',19980611);      %after passing rayleigh
    %    msg_awgn(j,:) = msg_ofdm(j,:);    %test, no noise

       %%%%%%%%%%%%%%%%%%test_no awgn%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    msg_nocp(j,:) = msg_ray(j,13:num_carrier+len_cp);    %wipe off cp
    %    msg_fft(j,:) = N*fft(msg_nocp(j,:));     %idft
       %%%%%%%%%%%%%%%DFT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       msg_nocp(j,:) = msg_awgn(j,13:num_carrier+len_cp);    %wipe off cp
       msg_fft(j,:) = N*fft(msg_nocp(j,:));     %idft

          %%%%%%%%%%%%%%%%%%理想信道估计%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       msg_nocp_ray(j,:) = msg_ray(j,13:num_carrier+len_cp);   %only rayleigh
       msg_nocp_est(j,:) = N*fft(msg_nocp_ray(j,:));

    end

    h_ray = msg_nocp_est./msg_prl;
    %%%%%%%%%%%%%%%%%%并串变换%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    msg_casd = reshape(transpose(msg_fft),1,length(msg_conv)/2);
    % H_RAY_casd = reshape(transpose(H_RAY),1,num_bit);
    h_ray_casd = reshape(transpose(h_ray),1,length(msg_conv)/2);

    %%%%%%%%%%%%%%%%%%%信道均衡(f)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    msg_casd_zf = msg_casd./h_ray_casd;
end

%%%%%%%%%%%%%%%%%%QPSK解调%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
msg_gry_demod = pskdemod(msg_casd_zf,M, pi/4);      %qpsk demodulate
[dummy, graydecode] = sort(grayencode); 
graydecode = graydecode - 1;
msg_demod = graydecode(msg_gry_demod+1)';
msg_demod_bi = de2bi(msg_demod,k)'; 
msg_demod_bi = msg_demod_bi(:);     %gray->normal

%%%%%%%%%%%%%%%%%%%解行列交织%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if more_num1 ~= 0
    msg_deintrlv_re = transpose(reshape(msg_demod_bi,row_intrlv*column_intrlv,group_num));      %row_intrlv*column_intrlvbits as a group again
    msg_deintrlv = zeros(group_num, row_intrlv*column_intrlv);
    for k = 1:group_num
        msg_deintrlv(k,:) = matdeintrlv(msg_deintrlv_re(k,:),row_intrlv, column_intrlv);   %relieve interleave
    end 
    msg_deintrlv_code = reshape(transpose(msg_deintrlv),1,size(msg_deintrlv,1)*size(msg_deintrlv,2));   %code after deinterleaving    
    deintrlv_len = length(msg_deintrlv_code);
    msg_deintrlv_code(deintrlv_len - rest_num1 + 1 : deintrlv_len) = [];     %eliminate 0
else
    msg_deintrlv_re = transpose(reshape(msg_demod_bi,row_intrlv*column_intrlv,group_num));      %row_intrlv*column_intrlvbits as a group again
    msg_deintrlv = zeros(group_num, row_intrlv*column_intrlv);
    for k = 1:group_num
        msg_deintrlv(k,:) = matdeintrlv(msg_deintrlv_re(k,:),row_intrlv, column_intrlv);   %relieve interleave
    end 
    msg_deintrlv_code = reshape(transpose(msg_deintrlv),1,size(msg_deintrlv,1)*size(msg_deintrlv,2));   %code after deinterleaving    
end

% %%%%%%%%%%%%%%%%%%解随机交织%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if more_num1 ~= 0
%     msg_deintrlv_re = transpose(reshape(msg_demod_bi,intrl_frm, intrl_rand_ts));
%     msg_deintrlv = zeros(intrl_rand_ts, intrl_frm);
%     for k = 1:intrl_rand_ts
%         msg_deintrlv(k, :) = randdeintrlv(msg_deintrlv_re(k, :), k);
%     end
%     msg_deintrlv_code = reshape(transpose(msg_deintrlv),1,size(msg_deintrlv,1)*size(msg_deintrlv,2));   %code after deinterleaving    
%     deintrl_len = length(msg_deintrlv_code);
%     msg_deintrlv_code(deintrl_len - rest_num1 + 1 : deintrl_len) = [];     %去掉0
% else
%     msg_deintrlv_re = transpose(reshape(msg_demod_bi,intrl_frm, intrl_rand_ts));
%     msg_deintrlv = zeros(intrl_rand_ts, intrl_frm);
%     for k = 1:intrl_rand_ts
%         msg_deintrlv(k, :) = randdeintrlv(msg_deintrlv_re(k, :), k);
%     end
%     msg_deintrlv_code = reshape(transpose(msg_deintrlv),1,size(msg_deintrlv,1)*size(msg_deintrlv,2));   %code after deinterleaving        
% end


% %%%%%%%%%%%%%%%%%%%解卷积码%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% msg_rx = vitdec(msg_deintrlv_code, trellis, tblen, 'cont', 'hard');     %'hard' judgment
% msg_rx(1:35) = [];  %delete traceback daley



%%%%%%%%%%%%%%%%%%%解Turbo码%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Turbo_decoder = comm.TurboDecoder(...
    'TrellisStructure', poly2trellis(4, [13 15 17], 13),...
    'InterleaverIndices', intrlvrIndices,...
    'NumIterations', 6);
msg_rx = step(Turbo_decoder, transpose(msg_deintrlv_code));     %Turbo decode
msg_rx = transpose(msg_rx);

%%%%%%%%%%%%%%%%%%%计算BER&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
% msg(num_bit-34: num_bit) = [];      %for convolutional code
msg_wrong = abs(msg_rx - msg);      %for turbo code
BER = sum(msg_wrong)/num_bit;


