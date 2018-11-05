% Regenerating data for Meng 
% Look at data from noBH volunteer acquisitions
%

look_at_2d_data = 1;
look_at_3d_data = 0;

if (look_at_2d_data == 1)
    which_vol = 3; which_case = 400; load_data = 1; report_roi = 1; save_gifs = 0; save_ft_uw_data = 0;
    Nx=256; Ny=256; Nz=4; Nt=4; Ne=2; Ncyc_gre=10; Ncyc_epi=30; Nti=8;
    mag_lim = [0,60]; phs_lim = [-pi,pi]; wave_lim = [-0.5,0.5]; mu_lim = [0,8]; loss_lim = [-4,4];
    md = inline('uint8((wandl(x,diff(a),mean(a))-a(1))*(b(2)-b(1))/(a(2)-a(1))+b(1));','x','a','b');
    
    if (which_case > 100), Nloops = Nz; else Nloops = 1; end
    for kloop = 1:Nloops
        switch which_vol
            case 1
                d0 = 'D:\Data\data-150306_noBH_vols\e3231_analy\';
                fb = 'e3231'; xs = 1:255; ys = 62:214;
            case 2
                d0 = 'D:\Data\data-150306_noBH_vols\e3232_analy\';
                fb = 'e3232'; xs = 1:255; ys = 53:205;
            case 3
                %d0 = 'D:\Data\data-150414_noBH_vols\e3316_analy\';
                d0 = '\\R5016755\OpusPresentations\Kevin\for_Meng\data-150414_noBH_vols\e3316_analy\';
                fb = 'e3316'; xs = 1:255; ys = 53:205;
            case 4
                d0 = 'D:\Data\data-150414_noBH_vols\e3317_analy\';
                fb = 'e3317'; xs = 1:255; ys = 53:205;
            case 5
                d0 = 'D:\Data\data-150415_noBH_vols\e3321_analy\';
                fb = 'e3321'; xs = 1:255; ys = 53:205;
        end
        
        switch which_case
            case 1  % BH GRE
                fin = [d0,fb,'_BH_GRE.mat']; fcimg = [d0,fb,'_BH_GRE.cimg'];
                fname_out = [d0,'gifs\',fb,'_BH_GRE.gif']; fname_out_mean_mag = [d0,'gifs\',fb,'_BH_GRE_mean_mag.gif'];
                mag_lim = [0,60]; bh_flag = 1;
            case 2  % BH EPI
                fin = [d0,fb,'_BH_EPI.mat']; fcimg = [d0,fb,'_BH_EPI.cimg'];
                fname_out = [d0,'gifs\',fb,'_BH_EPI.gif']; fname_out_mean_mag = [d0,'gifs\',fb,'_BH_EPI_mean_mag.gif'];
                mag_lim = [0,60]; bh_flag = 1;
            case {3,300}  % noBH GRE
                if (which_case < 100)
                    slice = 4; ss = num2str(slice); Ncyc = Ncyc_gre;
                    fin = [d0,fb,'_noBH_GRE_s',ss,'.mat']; fcimg = [d0,fb,'_noBH_GRE_s',ss,'.cimg'];
                    fname_out = [d0,'gifs\',fb,'_noBH_GRE_s',ss,'.gif']; fname_out_mean_mag = [d0,'gifs\',fb,'_noBH_GRE_mean_mag_s',ss,'.gif'];
                    fname_out_each_t = [d0,'gifs\',fb,'_noBH_GRE_eacht_s',ss,'.gif'];
                    fname_out_mean_mu = [d0,'gifs\',fb,'_noBH_GRE_mean_mu_s',ss,'.gif'];
                    fname_out_mean_loss = [d0,'gifs\',fb,'_noBH_GRE_mean_loss_s',ss,'.gif'];
                else
                    slice = kloop; ss = num2str(slice); Ncyc = Ncyc_gre;
                    fin = [d0,fb,'_noBH_GRE_s',ss,'.mat']; fcimg = [d0,fb,'_noBH_GRE_s',ss,'.cimg'];
                    fname_out = [d0,'gifs\',fb,'_noBH_GRE.gif']; fname_out_mean_mag = [d0,'gifs\',fb,'_noBH_GRE_mean_mag.gif'];
                    fname_out_each_t = [d0,'gifs\',fb,'_noBH_GRE_eacht'];
                    fname_out_mean_mu = [d0,'gifs\',fb,'_noBH_GRE_mean_mu.gif'];
                    fname_out_mean_loss = [d0,'gifs\',fb,'_noBH_GRE_mean_loss.gif'];
                    fout_mean_mag_data = [d0,fb,'_noBH_GRE_mean.mag'];
                    fout_ft_phs_data = [d0,fb,'_noBH_GRE_ft.phs'];
                end
                mag_lim = [0,60]; bh_flag = 0;
            case {4,400}  % noBH EPI
                if (which_case < 100)
                    slice = 4; ss = num2str(slice); Ncyc = Ncyc_epi;
                    fin = [d0,fb,'_noBH_EPI_s',ss,'.mat']; fcimg = [d0,fb,'_noBH_EPI_s',ss,'.cimg'];
                    fname_out = [d0,'gifs\',fb,'_noBH_EPI_s',ss,'.gif']; fname_out_mean_mag = [d0,'gifs\',fb,'_noBH_EPI_mean_mag_s',ss,'.gif'];
                    fname_out_each_t = [d0,'gifs\',fb,'_noBH_EPI_eacht_s',ss,'.gif'];
                    fname_out_mean_mu = [d0,'gifs\',fb,'_noBH_EPI_mean_mu_s',ss,'.gif'];
                    fname_out_mean_loss = [d0,'gifs\',fb,'_noBH_EPI_mean_loss_s',ss,'.gif'];
                else
                    slice = kloop; ss = num2str(slice); Ncyc = Ncyc_epi;
                    fin = [d0,fb,'_noBH_EPI_s',ss,'.mat']; fcimg = [d0,fb,'_noBH_EPI_s',ss,'.cimg'];
                    fname_out = [d0,'gifs\',fb,'_noBH_EPI.gif']; fname_out_mean_mag = [d0,'gifs\',fb,'_noBH_EPI_mean_mag.gif'];
                    fname_out_each_t = [d0,'gifs\',fb,'_noBH_EPI_eacht'];
                    fname_out_mean_mu = [d0,'gifs\',fb,'_noBH_EPI_mean_mu.gif'];
                    fname_out_mean_loss = [d0,'gifs\',fb,'_noBH_EPI_mean_loss.gif'];
                    fout_mean_mag_data = [d0,fb,'_noBH_EPI_mean.mag'];
                    fout_ft_phs_data = [d0,fb,'_noBH_EPI_ft.phs'];
                end
                mag_lim = [0,60]; bh_flag = 0;
        end
        
        if (bh_flag == 1)
            if (load_data == 1)
                C = single(load_nd_data(fcimg,[Ny,Nx,Nz,Nt],'float32','l',0,1));
                load(fin);
            end
            
            Nxs = length(xs); Nys = length(ys);
            CB = single(repmat(checkerboard(2,Ny/4,Nx/4),[1,1,Nz]));
            CB = permute(reshape(permute(CB(ys,xs,:),[2,3,1]),[Nz*Nxs,Nys]),[2,1])>0.5;
            
            M = zoom_array(permute(reshape(permute(abs(C(ys,xs,:,:)),[2,3,1,4]),[Nz*Nxs,Nys,Nt]),[2,1,3]),[1,1,2]);
            P = zoom_array(permute(reshape(permute(angle(C(ys,xs,:,:)),[2,3,1,4]),[Nz*Nxs,Nys,Nt]),[2,1,3]),[1,1,2]);
            W = permute(reshape(permute(mmdi_signal(ys,xs,:,:),[2,3,1,4]),[Nz*Nxs,Nys,8]),[2,1,3]);
            conf_mask = permute(reshape(permute(mmdi_lap_conf(ys,xs,:),[2,3,1]),[Nz*Nxs,Nys]),[2,1])>=0.95;
            %cmu = permute(reshape(permute(mmdi_cmu(ys,xs,:),[2,3,1]),[Nz*Nxs,Nys]),[2,1]);
            cmu = repmat(permute(reshape(permute(mmdi_cmu(ys,xs,:),[2,3,1]),[Nz*Nxs,Nys]),[2,1]).*(conf_mask | CB),[1,1,Nti]);
            M_mean = permute(reshape(permute(mean_mag(ys,xs,:),[2,3,1]),[Nz*Nxs,Nys]),[2,1]);
            
            if (report_roi == 1)
                fin_vals = zeros([Nz,6],'single');
                for k=1:Nz
                    indx = find(drawn_roi_mask(:,:,k)==0); tmp = mmdi_cmu(:,:,k); muv = tmp(indx);
                    fin_vals(k,:) = [mean(abs(muv)),std(abs(muv)), mean(real(muv)),std(real(muv)), mean(imag(muv)),std(imag(muv))];
                end
                fin_vals
            end
            
            if (save_gifs == 1)
                clear T;
                cmap_list = {gray(85),awave(85),aaasmo(85)};
                kr=1; kc=1; T{kr,kc}.data = M; T{kr,kc}.range = mag_lim; T{kr,kc}.cmap_indx = 1;
                kr=2; kc=1; T{kr,kc}.data = P; T{kr,kc}.range = phs_lim; T{kr,kc}.cmap_indx = 1;
                kr=3; kc=1; T{kr,kc}.data = W; T{kr,kc}.range = wave_lim; T{kr,kc}.cmap_indx = 2;
                kr=4; kc=1; T{kr,kc}.data = abs(cmu); T{kr,kc}.range = mu_lim; T{kr,kc}.cmap_indx = 3;
                kr=5; kc=1; T{kr,kc}.data = imag(cmu); T{kr,kc}.range = loss_lim; T{kr,kc}.cmap_indx = 2;
                multi_rgb_gif(cmap_list,T,fname_out,0.1,1);
                
                wl_data = md(M_mean,mag_lim,[1,256]); imwrite(wl_data,gray(256),fname_out_mean_mag,'GIF');
            end
        else
            Niter = Ncyc*Nt-Nt+1;
            if (load_data == 1)
                C = zoom_array(single(load_nd_data(fcimg,[Ny,Nx,Niter,Nt],'float32','l',0,1)),[1,1,1,2]);
                load(fin);
            end
            
            Nxs = length(xs); Nys = length(ys);
            C = reshape(permute(C(ys,xs,:,:),[1,2,4,3]),[Nys,Nxs,Nti*Niter]);
            CB = single(repmat(checkerboard(2,Ny/4,Nx/4),[1,1,Niter*Nti])); CB = CB(ys,xs,:)>0.5;
            CB2 = single(checkerboard(2,Ny/4,Nx/4)); CB2 = CB2(ys,xs,:)>0.5;
            
            M = abs(C); P = angle(C);
            try
                W = reshape(permute(mmdi_signal(ys,xs,:,:),[1,2,4,3]),[Nys,Nxs,Nti*Niter]);
            catch
                mmdi_signal = load_nd_data([fcimg(1:(end-5)),'.mmdi_signal'],[Ny,Nx,Niter,Nti],'float32','l',0,0);
                W = reshape(permute(mmdi_signal(ys,xs,:,:),[1,2,4,3]),[Nys,Nxs,Nti*Niter]);
            end
            conf_mask = zoom_array(mmdi_lap_conf(ys,xs,:),[1,1,Nti])>=0.95;
            %cmu = permute(reshape(permute(mmdi_cmu(ys,xs,:),[2,3,1]),[Nz*Nxs,Nys]),[2,1]);
            cmu = zoom_array(mmdi_cmu(ys,xs,:),[1,1,Nti]).*(conf_mask | CB);
            M_mean = mean_mag(ys,xs);
            mean_cmu = mean(mmdi_cmu(ys,xs,:),3).*((mean_lap_conf(ys,xs)>=0.95) | CB2);
            if (save_ft_uw_data == 1)
                tmp = cat(3,mmdi_uw(:,:,:,1),squeeze(mmdi_uw(:,:,end,2:end)));
                ft_uw = ifft(tmp,[],3); h1_uw = ft_uw(:,:,Ncyc+1);
                ts = repmat(reshape(0:(Nt-1),[1,1,Nt]),[Ny,Nx,1]);
                h1t_uw = 2*repmat(abs(h1_uw),[1,1,Nt]).*sin(repmat(angle(h1_uw),[1,1,Nt])-2*pi*ts/Nt);
            end
            
            if (report_roi == 1)
                if (which_case < 100)
                    fin_vals = zeros([Niter,6],'single');
                    indx = find(drawn_roi_mask==0);
                    for k=1:Niter
                        tmp = mmdi_cmu(:,:,k); muv = tmp(indx);
                        fin_vals(k,:) = [mean(abs(muv)),std(abs(muv)), mean(real(muv)),std(real(muv)), mean(imag(muv)),std(imag(muv))];
                    end
                    %fin_vals;
                else
                    if (slice == 1), fin_vals = zeros([Niter,6,Nz],'single'); end
                    indx = find(drawn_roi_mask==0);
                    for k=1:Niter
                        tmp = mmdi_cmu(:,:,k); muv = tmp(indx);
                        fin_vals(k,:,slice) = [mean(abs(muv)),std(abs(muv)), mean(real(muv)),std(real(muv)), mean(imag(muv)),std(imag(muv))];
                    end
                    %fin_vals;
                    if (slice == Nz)
                        for kplot = 1:2
                            if (kplot == 1)
                                ystr = 'Stiffness (kPa)'; mindx=1; sdindx=2;
                            elseif (kplot == 2)
                                ystr = 'Loss Modulus (kPa)'; mindx=5; sdindx=6;
                            end
                            x_lim = ceil(Niter/10)*10;
                            hfig = figure; colordef_kjg(hfig,'word2'); set(hfig,'Position',[20,80,1800,900]);
                            for ksl = 1:Nz
                                mvals = fin_vals(:,mindx,ksl); sdvals = fin_vals(:,sdindx,ksl);
                                m_of_m = mean(mvals); sd_of_m = std(mvals);
                                subplot(2,2,ksl); errorbar((1:Niter)',mvals,sdvals); grid on; axis([0,x_lim,0,4]);
                                ylabel(ystr); xlabel('Trial #'); title(sprintf('Slice %d ROIs: %0.3f +- %0.3f',ksl,m_of_m,sd_of_m));
                            end
                        end
                    end
                end
            end
            
            if (save_gifs == 1)
                if (which_case < 100)
                    clear T;
                    cmap_list = {gray(85),awave(85),aaasmo(85)};
                    kr=1; kc=1; T{kr,kc}.data = M; T{kr,kc}.range = mag_lim; T{kr,kc}.cmap_indx = 1;
                    kr=2; kc=1; T{kr,kc}.data = P; T{kr,kc}.range = phs_lim; T{kr,kc}.cmap_indx = 1;
                    kr=3; kc=1; T{kr,kc}.data = W; T{kr,kc}.range = wave_lim; T{kr,kc}.cmap_indx = 2;
                    kr=4; kc=1; T{kr,kc}.data = abs(cmu); T{kr,kc}.range = mu_lim; T{kr,kc}.cmap_indx = 3;
                    kr=5; kc=1; T{kr,kc}.data = imag(cmu); T{kr,kc}.range = loss_lim; T{kr,kc}.cmap_indx = 2;
                    multi_rgb_gif(cmap_list,T,fname_out,0.1,1);
                    
                    wl_data = md(M_mean,mag_lim,[1,256]); imwrite(wl_data,gray(256),fname_out_mean_mag,'GIF');
                    wl_data = md(mean(abs(cmu),3),mu_lim,[1,256]); imwrite(wl_data,aaasmo(256),fname_out_mean_mu,'GIF');
                    wl_data = md(mean(imag(cmu),3),loss_lim,[1,256]); imwrite(wl_data,awave(256),fname_out_mean_loss,'GIF');
                    
                    clear T;
                    cmap_list = {gray(85),awave(85),aaasmo(85)};
                    kr=1; kc=1; T{kr,kc}.data = cat(2,M(:,:,1:(Nti*Nt):end),M(:,:,3:(Nti*Nt):end),M(:,:,5:(Nti*Nt):end),M(:,:,7:(Nti*Nt):end)); T{kr,kc}.range = mag_lim; T{kr,kc}.cmap_indx = 1;
                    kr=2; kc=1; T{kr,kc}.data = cat(2,P(:,:,1:(Nti*Nt):end),P(:,:,3:(Nti*Nt):end),P(:,:,5:(Nti*Nt):end),P(:,:,7:(Nti*Nt):end)); T{kr,kc}.range = phs_lim; T{kr,kc}.cmap_indx = 1;
                    kr=3; kc=1; T{kr,kc}.data = cat(2,W(:,:,1:(Nti*Nt):end),W(:,:,3:(Nti*Nt):end),W(:,:,5:(Nti*Nt):end),W(:,:,7:(Nti*Nt):end)); T{kr,kc}.range = wave_lim; T{kr,kc}.cmap_indx = 2;
                    kr=4; kc=1; T{kr,kc}.data = cat(2,abs(cmu(:,:,1:(Nti*Nt):end)),abs(cmu(:,:,3:(Nti*Nt):end)),abs(cmu(:,:,5:(Nti*Nt):end)),abs(cmu(:,:,7:(Nti*Nt):end))); T{kr,kc}.range = mu_lim; T{kr,kc}.cmap_indx = 3;
                    kr=5; kc=1; T{kr,kc}.data = cat(2,imag(cmu(:,:,1:(Nti*Nt):end)),imag(cmu(:,:,3:(Nti*Nt):end)),imag(cmu(:,:,5:(Nti*Nt):end)),imag(cmu(:,:,7:(Nti*Nt):end))); T{kr,kc}.range = loss_lim; T{kr,kc}.cmap_indx = 2;
                    multi_rgb_gif(cmap_list,T,fname_out_each_t,0.1,1);
                else
                    if (slice == 1), clear T; end
                    cmap_list = {gray(85),awave(85),aaasmo(85)};
                    kr=1; kc=slice; T{kr,kc}.data = M; T{kr,kc}.range = mag_lim; T{kr,kc}.cmap_indx = 1;
                    kr=2; kc=slice; T{kr,kc}.data = P; T{kr,kc}.range = phs_lim; T{kr,kc}.cmap_indx = 1;
                    kr=3; kc=slice; T{kr,kc}.data = W; T{kr,kc}.range = wave_lim; T{kr,kc}.cmap_indx = 2;
                    kr=4; kc=slice; T{kr,kc}.data = abs(cmu); T{kr,kc}.range = mu_lim; T{kr,kc}.cmap_indx = 3;
                    kr=5; kc=slice; T{kr,kc}.data = imag(cmu); T{kr,kc}.range = loss_lim; T{kr,kc}.cmap_indx = 2;
                    if (slice == Nz), multi_rgb_gif(cmap_list,T,fname_out,0.1,1); end
                    
                    if (slice == 1)
                        %Mv = M_mean; Amuv = mean(abs(cmu),3); Imuv = mean(imag(cmu),3);
                        Mv = M_mean; Amuv = abs(mean_cmu); Imuv = imag(mean_cmu);
                    else
                        %Mv = cat(2,Mv,M_mean); Amuv = cat(2,Amuv,mean(abs(cmu),3)); Imuv = cat(2,Imuv,mean(imag(cmu),3));
                        Mv = cat(2,Mv,M_mean); Amuv = cat(2,Amuv,abs(mean_cmu)); Imuv = cat(2,Imuv,imag(mean_cmu));
                    end
                    if (slice == Nz)
                        wl_data = md(Mv,mag_lim,[0,255]); imwrite(wl_data,gray(256),fname_out_mean_mag,'GIF');
                        wl_data = md(Amuv,mu_lim,[0,255]); imwrite(wl_data,aaasmo(256),fname_out_mean_mu,'GIF');
                        wl_data = md(Imuv,loss_lim,[0,255]); imwrite(wl_data,awave(256),fname_out_mean_loss,'GIF');
                    end
                    
                    clear T2;
                    cmap_list = {gray(85),awave(85),aaasmo(85)};
                    kr=1; kc=1; T2{kr,kc}.data = cat(2,M(:,:,1:(Nti*Nt):end),M(:,:,3:(Nti*Nt):end),M(:,:,5:(Nti*Nt):end),M(:,:,7:(Nti*Nt):end)); T2{kr,kc}.range = mag_lim; T2{kr,kc}.cmap_indx = 1;
                    kr=2; kc=1; T2{kr,kc}.data = cat(2,P(:,:,1:(Nti*Nt):end),P(:,:,3:(Nti*Nt):end),P(:,:,5:(Nti*Nt):end),P(:,:,7:(Nti*Nt):end)); T2{kr,kc}.range = phs_lim; T2{kr,kc}.cmap_indx = 1;
                    kr=3; kc=1; T2{kr,kc}.data = cat(2,W(:,:,1:(Nti*Nt):end),W(:,:,3:(Nti*Nt):end),W(:,:,5:(Nti*Nt):end),W(:,:,7:(Nti*Nt):end)); T2{kr,kc}.range = wave_lim; T2{kr,kc}.cmap_indx = 2;
                    kr=4; kc=1; T2{kr,kc}.data = cat(2,abs(cmu(:,:,1:(Nti*Nt):end)),abs(cmu(:,:,3:(Nti*Nt):end)),abs(cmu(:,:,5:(Nti*Nt):end)),abs(cmu(:,:,7:(Nti*Nt):end))); T2{kr,kc}.range = mu_lim; T2{kr,kc}.cmap_indx = 3;
                    kr=5; kc=1; T2{kr,kc}.data = cat(2,imag(cmu(:,:,1:(Nti*Nt):end)),imag(cmu(:,:,3:(Nti*Nt):end)),imag(cmu(:,:,5:(Nti*Nt):end)),imag(cmu(:,:,7:(Nti*Nt):end))); T2{kr,kc}.range = loss_lim; T2{kr,kc}.cmap_indx = 2;
                    multi_rgb_gif(cmap_list,T2,[fname_out_each_t,'_s',num2str(slice),'.gif'],0.1,1);
                end
            end
            
            if (save_ft_uw_data == 1)
                if (which_case >= 100)
                    if (slice == 1), T_ft = zeros([Ny,Nx,Nz,Nt],'single'); M_means = zeros([Ny,Nx,Nz],'single'); end
                    T_ft(:,:,slice,:) = permute(h1t_uw,[1,2,4,3]); M_means(:,:,slice) = mean_mag;
                    if (slice == Nz)
                        save_nd_data(T_ft,fout_ft_phs_data,'float32','l',0,0);
                        save_nd_data(M_means,fout_mean_mag_data,'float32','l',0,0);
                    end
                end
            end
        end  % if (bh_flag == 1)
    end  % for kloop = 1:Nloops
end

if (look_at_3d_data == 1)
    organize_mre_data = 1;
    make_ppt_file = 1;
    
    proc_list = []; dproc = []; dalt = [];
    if (1==0)
        %dbase = 'D:\Data\data-150414_noBH_vols\e3316_s5_21757_15_04_14\';
        dbase = 'D:\Data\data-150414_noBH_vols\e3316_s6_11958_15_04_14\';
        dout = 'D:\Data\data-150414_noBH_vols\e3316_summary\';
        dexam = dbase;
        proc_list = {dbase};
        %hfile = [dbase,'E3316S5I1.MR']; bser = 5; MRE_data_series = bser; % 60-Hz
        hfile = [dbase,'E3316S6I1.MR']; bser = 6; MRE_data_series = bser; % 60-Hz
        size_info = [256,256,32,3,2]; Nti=9;
    elseif (1==0)
        %dbase = 'D:\Data\data-150414_noBH_vols\e3317_s5_3480_15_04_14\';
        dbase = 'D:\Data\data-150414_noBH_vols\e3317_s6_19317_15_04_14\';
        dout = 'D:\Data\data-150414_noBH_vols\e3317_summary\';
        dexam = dbase;
        proc_list = {dbase};
        %hfile = [dbase,'E3317S5I1.MR']; bser = 5; MRE_data_series = bser; % 60-Hz
        hfile = [dbase,'E3317S6I1.MR']; bser = 6; MRE_data_series = bser; % 60-Hz
        size_info = [256,256,32,3,2]; Nti=9;
    elseif (1==1)
        %dbase = 'D:\Data\data-150415_noBH_vols\e3321_s5_1386_15_04_15\';
        dbase = 'D:\Data\data-150415_noBH_vols\e3321_s6_20062_15_04_15\';
        dout = 'D:\Data\data-150415_noBH_vols\e3321_summary\';
        dexam = dbase;
        proc_list = {dbase};
        %hfile = [dbase,'E3321S5I1.MR']; bser = 5; MRE_data_series = bser; % 60-Hz
        hfile = [dbase,'E3321S6I1.MR']; bser = 6; MRE_data_series = bser; % 60-Hz
        size_info = [256,256,32,3,2]; Nti=9;
    end
    
    if (organize_mre_data == 1)
        %ser_list = 1:13;
        ser_list = 1:15;
        for kproc = 1:length(proc_list)
            dproc = []; cur_proc = dbase;
            
            for kser = ser_list
                if (isempty(dproc))
                    ser_str = num2str(bser*100+kser);
                    %dser = [dexam,'s',ser_str,'\']; D = dir([dser,'i*.dcm']); Nfiles = length(D);
                    
                    %Dser = dir([dexam,'*_',ser_str]); dser = [dexam,Dser(1).name,'\'];
                    %D = dir([dser,'IM*.dcm']); Nfiles = length(D);
                    
                    dser = dexam; D = dir([dser,'epidcm_s',num2str(kser),'i*.proc']); Nfiles = length(D);
                end
                for k = 1:Nfiles
                    fname = [dser,D(k).name];
                    if (k==1)
                        H = dicominfo(fname); ser_str = num2str(H.SeriesNumber);
                        I = dicomread(fname); Ny=size(I,1); Nx=size(I,2); Nz=size(I,3);
                        S = zeros([Ny,Nx,Nz,Nfiles],'int16');
                    end
                    S(:,:,:,k) = dicomread(fname);
                end
                Tstr = sprintf('%dx',size(S)); Tstr = ['_',Tstr(1:(end-1)),'.img'];
                save_nd_data(S,[dout,'S',ser_str,Tstr],'int16','l',0,0);
            end
        end
    end
    
    if (make_ppt_file == 1)
        load_files = 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Location and name of files
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        gif_dir = 'D:\Data\gif_dir\';  % output GIF folder
        
        cbar_bmp_dir = 'D:\Data\gif_dir\';
        black_box_fname = [cbar_bmp_dir,'black_box.bmp'];
        gray_cbar_fname = [cbar_bmp_dir,'gray_cbar.bmp'];
        awave_cbar_fname = [cbar_bmp_dir,'awave_cbar.bmp'];
        aaasmo_cbar_fname = [cbar_bmp_dir,'aaasmo_cbar.bmp'];
        %pptx_fname_out = [dbase,'pancreas_mre_summary_autogen.pptx'];
        pptx_fname_out = [dout,'mre_summary_autogen.pptx'];
        Nx=size_info(1); Ny=size_info(2); Nz=size_info(3); Nt = size_info(4);
        if (~isempty(Nti) && (Nti>Nt)), Nti = ceil(Nti/Nt)*Nt; else Nti = Nt; end
        
        %%%%%%%%%%%%%%%%%%%%%
        % Load in data
        %%%%%%%%%%%%%%%%%%%%%
        if (load_files == 1)
            if (~isempty(dalt))
                use_mean_mag_flag = 0;
                din = dalt;
                S1 = permute(single(load_nd_data([din,'epimre_x.phs'],[Ny,Nx,Nz,Nt],'float32','l',0,0)),[1,2,4,3]);
                S2 = permute(single(load_nd_data([din,'epimre_y.phs'],[Ny,Nx,Nz,Nt],'float32','l',0,0)),[1,2,4,3]);
                S3 = permute(single(load_nd_data([din,'epimre_z.phs'],[Ny,Nx,Nz,Nt],'float32','l',0,0)),[1,2,4,3]);
                tmp_ft = zeros([Ny,Nx,Nt,Nz],'single');
                tmp = single(load_nd_data([din,'epimre_spua.cx'],[Ny,Nx,1,Nz],'float32','l',0,1));
                tmp_ft(:,:,2,:) = tmp; tmp_ft(:,:,end,:)=conj(tmp); S4 = real(fft(tmp_ft,[],3));
                tmp = single(load_nd_data([din,'epimre_spua.cy'],[Ny,Nx,1,Nz],'float32','l',0,1));
                tmp_ft(:,:,2,:) = tmp; tmp_ft(:,:,end,:)=conj(tmp); S5 = real(fft(tmp_ft,[],3));
                tmp = single(load_nd_data([din,'epimre_spua.cz'],[Ny,Nx,1,Nz],'float32','l',0,1));
                tmp_ft(:,:,2,:) = tmp; tmp_ft(:,:,end,:)=conj(tmp); S6 = real(fft(tmp_ft,[],3));
                tmp = single(load_nd_data([din,'epimre_spua.div'],[Ny,Nx,1,Nz],'float32','l',0,1));
                tmp_ft(:,:,2,:) = tmp; tmp_ft(:,:,end,:)=conj(tmp); S7 = real(fft(tmp_ft,[],3));
                clear tmp_ft
                max_val = max([max(S4(:)),max(S5(:)),max(S6(:)),max(S7(:))]);
                S4=S4*32000/max_val; S5=S5*32000/max_val; S6=S6*32000/max_val; S7=S7*32000/max_val;
                %tmp = single(load_nd_data([din,'epimre_spua_3x3x3di_mf.cmu'],[Ny,Nx,1,Nz],'float32','l',0,1));
                tmp = single(load_nd_data([din,'epimre_spua_5x5x3di_mf.cmu'],[Ny,Nx,1,Nz],'float32','l',0,1));
                S8 = abs(tmp); S9 = real(tmp); S10 = imag(tmp);
                if (use_mean_mag_flag == 1)
                    S14 = single(load_nd_data([din,'epimre_prezf_mean.mag'],[Ny,Nx,1,Nz],'float32','l',0,0));
                else
                    S14 = permute(single(load_nd_data([din,'epimre_prezf_x.mag'],[Ny,Nx,Nz,Nt],'float32','l',0,0)),[1,2,4,3]);
                end
                S14 = S14*32000/max(S14(:));
                S15 = ones(size(S8),'single');
                din = dout;
            else
                use_mean_mag_flag = 0; bser = MRE_data_series;
                %din = dbase; %Nx=128; Ny=128; Nz=48; Nt=8;
                din = dout;
                % 1,2,3 = Px, Py, Pz
                % 4,5,6,7 = Cx, Cy, Cz, Div
                % 8,9,10 = G, Gr, Gi
                % 14,15 = mag, lfe conf map
                for k = [1:10,14,15]
                    switch k
                        case {1,2,3}
                            fin_size = [Ny,Nx,Nt,Nz]; sf=1000;
                        case {4,5,6,7}
                            fin_size = [Ny,Nx,Nt,Nz]; sf=1;
                        case {8,9,10,15}
                            fin_size = [Ny,Nx,1,Nz]; sf=1000;
                        case {14}
                            fin_size = [Ny,Nx,1,Nz]; sf=1;
                    end
                    eval(['S',num2str(k),' = [];']);
                    D = dir([din,'S',num2str(bser*100+k),'*.img']);
                    if (~isempty(D))
                        T = single(load_nd_data([din,D(1).name],fin_size,'int16','l',0,0))/sf;
                        if (~isempty(Nti) && (Nti > Nt))
                            switch k
                                case {1,2,3}
                                    T = zoom_array(T,[1,1,Nti/Nt,1]);
                                case {4,5,6,7}
                                    tmp = ifft(T,[],3); T=zeros([Ny,Nx,Nti,Nz],'single');
                                    T(:,:,2,:) = tmp(:,:,2,:); T(:,:,end,:) = conj(tmp(:,:,2,:)); T = real(fft(T,[],3));
                                case {8,9,10,15}
                                    % Nothing
                                case {14}
                                    % Nothing here
                            end
                        end
                        eval(['S',num2str(k),' = T;']);
                        if (k==14), S14 = repmat(S14,[1,1,Nti,1]); end
                    end
                end
                if (use_mean_mag_flag == 0)
                    D = dir([din,'S',num2str(bser*100+1),'*.img']);
                    T = load_nd_data([din,D(1).name],[Ny,Nx,Nt,Nz,2],'int16','l',0,0);
                    S14 = single(T(:,:,:,:,2));
                    if (~isempty(Nti) && (Nti > Nt)), S14 = zoom_array(S14,[1,1,Nti/Nt,1]); end
                end
                if ((isempty(S15)) && (~isempty(S8))), S15 = ones(size(S8),'single'); end
            end
            
            %D = dir([din,'resampled_T1.vol']);
            try
                D = dir([din,'resampled_S',num2str(T1_series),'.vol']);
            catch
                D = [];
            end
            if (isempty(D))
                %T1vol = zeros([Ny,Nx,Nt,Nz],'single');
                T1vol = [];
                %T1vol = S7;
            else
                T1vol = single(load_nd_data([din,D(1).name],[Ny,Nx,1,Nz],'float32','l',0,0));
                T1vol = repmat(T1vol,[1,1,Nti,1]);
            end
            
            if (~isempty(hfile))
                H = dicominfo(hfile);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Get information from the DICOM header to put into Powerpoint file.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (isfield(H,'PatientName'))
                    PatientName = H.PatientName;
                    if (isfield(PatientName,'FamilyName'))
                        FamilyName = PatientName.FamilyName;
                        if (isempty(FamilyName)), FamilyName = 'X'; end
                    else
                        FamilyName = 'X';
                    end
                else
                    FamilyName = 'X';
                end
                
                if (isfield(H,'PatientID'))
                    PatientID = H.PatientID;
                    if (isempty(PatientID)), PatientID = '000'; end
                    if (length(PatientID)<3), PatientID = '000'; end
                else
                    PatientID = '000';
                end
                
                PatientCode = [FamilyName(1),PatientID((end-2):end)];  % e.g., 'H532'
                if (isfield(H,'PatientSex')), PatientSex = H.PatientSex; else PatientSex = 'M/F?'; end  % e.g., 'M'
                if (isfield(H,'PatientAge')), PatientAge = H.PatientAge; else PatientAge = '???Y'; end  % e.g., '053Y'
                if (isfield(H,'StudyDate')), StudyDate = H.StudyDate; else StudyDate = '20XXXXXX'; end  % e.g., '20130916'
                if (isfield(H,'StudyID')), StudyID = H.StudyID; else StudyID = 'XXXXX'; end  % e.g., '10704'
                if (isfield(H,'SeriesNumber')), Series = H.SeriesNumber; else Series = '0'; end  % e.g., 5
            else
                FamilyName = 'X'; PatientID = '000'; PatientCode = [FamilyName(1),PatientID((end-2):end)];
                PatientSex = 'M/F?'; PatientAge = '???Y'; StudyDate = '20XXXXXX'; StudyID = 'XXXXX'; Series = '0';
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create animated GIFs
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Set limits for windowing and leveling the data for display purposes.
        mag_lim = [0,32000]; phs_lim = [-pi,pi]; curl_lim = [-32000,32000];
        Gr_lim = [0,8]; Gi_lim = [-4,4]; T1_lim = [0,8000];
        cmap_list = {gray(64),awave(64),aaasmo(64)};
        add_slice_information = 1;
        mask_data_flag = 1; conf_thresh = 0;
        
        if (mask_data_flag == 1)
            conf_mask = S15>=conf_thresh;
            CB = repmat(checkerboard(2,Ny/4,Nx/4),[1,1,1,Nz])>0.5;
        else
            conf_mask = ones(size(S14))==1; CB = ones(size(S14))==1;
        end
        S8 = S8.*(conf_mask|CB); S9 = S9.*(conf_mask|CB); S10 = S10.*(conf_mask|CB);
        
        ani_list = [1,(1:Nz)+1];
        fname_cell = cell(length(ani_list),1);
        for which_ani = ani_list
            switch which_ani
                case 1
                    zrng = 1:Nz; fname_out = [gif_dir,'full_brain_volume.gif'];
                otherwise
                    zrng = which_ani-1; fname_out = [gif_dir,sprintf('brain_slice_%0.2d.gif',zrng)];
            end
            fname_cell{which_ani,1} = fname_out;
            % Save the whole volume as an animated GIF
            clear S
            %S = struct('data',cell(3,4),'range',cell(3,4),'cmap_indx',cell(3,4));
            S = cell(3,4);
            kr=1; kc=1; S{kr,kc}.data = S14(:,:,:,zrng); S{kr,kc}.range = mag_lim; S{kr,kc}.cmap_indx = 1;
            kr=1; kc=2; S{kr,kc}.data = S1(:,:,:,zrng); S{kr,kc}.range = phs_lim; S{kr,kc}.cmap_indx = 1;
            kr=1; kc=3; S{kr,kc}.data = S2(:,:,:,zrng); S{kr,kc}.range = phs_lim; S{kr,kc}.cmap_indx = 1;
            kr=1; kc=4; S{kr,kc}.data = S3(:,:,:,zrng); S{kr,kc}.range = phs_lim; S{kr,kc}.cmap_indx = 1;
            if (isempty(T1vol))
                kr=2; kc=1; S{kr,kc}.data = S7(:,:,:,zrng); S{kr,kc}.range = curl_lim; S{kr,kc}.cmap_indx = 2;
            else
                kr=2; kc=1; S{kr,kc}.data = T1vol(:,:,:,zrng); S{kr,kc}.range = T1_lim; S{kr,kc}.cmap_indx = 1;
            end
            kr=2; kc=2; S{kr,kc}.data = S4(:,:,:,zrng); S{kr,kc}.range = curl_lim; S{kr,kc}.cmap_indx = 2;
            kr=2; kc=3; S{kr,kc}.data = S5(:,:,:,zrng); S{kr,kc}.range = curl_lim; S{kr,kc}.cmap_indx = 2;
            kr=2; kc=4; S{kr,kc}.data = S6(:,:,:,zrng); S{kr,kc}.range = curl_lim; S{kr,kc}.cmap_indx = 2;
            if (add_slice_information == 1)
                w = 32000;  % value to get white numbers using mag_lim window-and-level and gray colormap
                T = 0*S14(:,:,:,zrng); z = 4; Np = 5*z; b = repmat(single(0),[2*Np,z]); cnt = 0;
                for kz = zrng
                    cnt = cnt + 1;
                    digs = [floor(kz/10),mod(kz,10)];
                    digstr = cat(2,b,block_number(digs(1),z),b,block_number(digs(2),z));
                    size_digs = size(digstr);
                    for kt = 1:Nti
                        T((end-z-(size_digs(1)-1)):(end-z),1:size_digs(2),kt,cnt) = digstr*w;
                    end
                end
            else
                T = 0*S14(:,:,:,zrng);
            end
            kr=3; kc=1; S{kr,kc}.data = T; S{kr,kc}.range = mag_lim; S{kr,kc}.cmap_indx = 1;
            kr=3; kc=2; S{kr,kc}.data = repmat(S8(:,:,:,zrng),[1,1,Nti,1]); S{kr,kc}.range = Gr_lim; S{kr,kc}.cmap_indx = 3;
            kr=3; kc=3; S{kr,kc}.data = repmat(S9(:,:,:,zrng),[1,1,Nti,1]); S{kr,kc}.range = Gr_lim; S{kr,kc}.cmap_indx = 3;
            kr=3; kc=4; S{kr,kc}.data = repmat(S10(:,:,:,zrng),[1,1,Nti,1]); S{kr,kc}.range = Gi_lim; S{kr,kc}.cmap_indx = 2;
            
            multi_rgb_gif(cmap_list,S,fname_out,0.1,1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Create Powerpoint presentation with animations in it.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ws_flag = 0; % widescreen presentation flag
        if (ws_flag == 1)
            ppt_dim = [13.333, 7.5];  % Powerpoint slide size in inches
        else
            ppt_dim = [10, 7.5];
        end
        
        % Start new presentation
        isOpen  = exportToPPTX();
        if ~isempty(isOpen),
            % If PowerPoint already started, then close first and then open a new one
            exportToPPTX('close');
        end
        
        exportToPPTX('new','Dimensions',ppt_dim, ...
            'Title','MR Elastography Summary', ...
            'Author',' ', ...
            'Subject','MRE Summary', ...
            'Comments','Summary of MRE results.');
        exportToPPTX('addslide');
        exportToPPTX('addtext',[PatientCode,': ',PatientAge,' ',PatientSex],'Position',[0, ppt_dim(2)/2-3, ppt_dim(1), 3],'Color',[0,0,0],'BackgroundColor',[1,1,1],...
            'FontSize',40,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
        %exportToPPTX('addtext',['Date: ',StudyDate,', Exam ',StudyID],'Position',[0, ppt_dim(2)/2, ppt_dim(1), 3],'Color',[0,0,0],'BackgroundColor',[1,1,1],...
        %    'FontSize',40,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
        exportToPPTX('addtext',['Date: ',StudyDate,', Exam ',StudyID,', S',num2str(Series)],'Position',[0, ppt_dim(2)/2, ppt_dim(1), 3],'Color',[0,0,0],'BackgroundColor',[1,1,1],...
            'FontSize',40,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
        
        for ksl = 1:size(fname_cell,1)
            fname_in = fname_cell{ksl,1}; make_cmap_flag = 1;
            w = 10; h = 7.5; tbh = 0.25;
            exportToPPTX('addslide');
            exportToPPTX('addpicture',fname_in,'Position',[0, 0, w, h]);
            
            exportToPPTX('addtext','Magnitude','Position',[0, 0, w/4, tbh],'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                'FontSize',16,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
            exportToPPTX('addtext','X Phase','Position',[w/4, 0, w/4, tbh],'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                'FontSize',16,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
            exportToPPTX('addtext','Y Phase','Position',[2*w/4, 0, w/4, tbh],'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                'FontSize',16,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
            exportToPPTX('addtext','Z Phase','Position',[3*w/4, 0, w/4, tbh],'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                'FontSize',16,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
            
            if (isempty(T1vol))
                exportToPPTX('addtext','Divergence','Position',[0, h/3, w/4, tbh],'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                    'FontSize',16,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
            else
                exportToPPTX('addtext','3D MP-RAGE','Position',[0, h/3, w/4, tbh],'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                    'FontSize',16,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
            end
            exportToPPTX('addtext','X Curl','Position',[w/4, h/3, w/4, tbh],'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                'FontSize',16,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
            exportToPPTX('addtext','Y Curl','Position',[2*w/4, h/3, w/4, tbh],'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                'FontSize',16,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
            exportToPPTX('addtext','Z Curl','Position',[3*w/4, h/3, w/4, tbh],'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                'FontSize',16,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
            
            exportToPPTX('addtext','Stiffness','Position',[w/4, 2*h/3, w/4, tbh],'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                'FontSize',16,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
            exportToPPTX('addtext','Storage Modulus','Position',[2*w/4, 2*h/3, w/4, tbh],'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                'FontSize',16,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
            exportToPPTX('addtext','Loss Modulus','Position',[3*w/4, 2*h/3, w/4, tbh],'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                'FontSize',16,'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle');
            
            if (make_cmap_flag == 1)
                cmap_font_size = 8; tbh = 0.14;
                bb_vec = [w/4, 2*h/3+0.66*h/3, 0.1875*w/4, 0.33*h/3];
                map_pos = bb_vec(1:2)+[0,0.05]; map_size = [bb_vec(3)/3, bb_vec(4)-0.05];  % [H,V], [W,H]
                top_pos_vec = [bb_vec(1)+map_size(1), bb_vec(2), 0.66*bb_vec(3), tbh];
                bot_pos_vec = [bb_vec(1)+map_size(1), bb_vec(2)+bb_vec(4)-tbh, 0.66*bb_vec(3), tbh];
                
                exportToPPTX('addpicture',black_box_fname,'Position',bb_vec);
                exportToPPTX('addpicture',aaasmo_cbar_fname,'Position',[map_pos, map_size]);
                exportToPPTX('addtext',sprintf('%d',Gr_lim(2)),'Position',top_pos_vec,'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                    'FontSize',cmap_font_size,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
                exportToPPTX('addtext',sprintf('%d',Gr_lim(1)),'Position',bot_pos_vec,'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                    'FontSize',cmap_font_size,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
                
                bb_vec = bb_vec + [w/4, 0, 0, 0];
                map_pos = map_pos + [w/4, 0];  % [H,V], [W,H]
                top_pos_vec = top_pos_vec + [w/4, 0, 0, 0];
                bot_pos_vec = bot_pos_vec + [w/4, 0, 0, 0];
                
                exportToPPTX('addpicture',black_box_fname,'Position',bb_vec);
                exportToPPTX('addpicture',aaasmo_cbar_fname,'Position',[map_pos, map_size]);
                exportToPPTX('addtext',sprintf('%d',Gr_lim(2)),'Position',top_pos_vec,'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                    'FontSize',cmap_font_size,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
                exportToPPTX('addtext',sprintf('%d',Gr_lim(1)),'Position',bot_pos_vec,'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                    'FontSize',cmap_font_size,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
                
                
                bb_vec = bb_vec + [w/4, 0, 0, 0];
                map_pos = map_pos + [w/4, 0];  % [H,V], [W,H]
                top_pos_vec = top_pos_vec + [w/4, 0, 0, 0];
                bot_pos_vec = bot_pos_vec + [w/4, 0, 0, 0];
                
                exportToPPTX('addpicture',black_box_fname,'Position',bb_vec);
                exportToPPTX('addpicture',awave_cbar_fname,'Position',[map_pos, map_size]);
                exportToPPTX('addtext',sprintf('%d',Gi_lim(2)),'Position',top_pos_vec,'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                    'FontSize',cmap_font_size,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
                exportToPPTX('addtext',sprintf('%d',Gi_lim(1)),'Position',bot_pos_vec,'Color',[1,1,1],'BackgroundColor',[0,0,0],...
                    'FontSize',cmap_font_size,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','middle');
            end
        end
        
        % Check current presentation
        fileStats = exportToPPTX('query');
        
        if ~isempty(fileStats),
            fprintf('Presentation size: %f x %f\n',fileStats.dimensions);
            fprintf('Number of slides: %d\n',fileStats.numSlides);
        end
        
        % Save presentation -- overwrite file if it already exists
        % Filename automatically checked for proper extension
        newFile = exportToPPTX('save',pptx_fname_out);
        % Close presentation (and clear all temporary files)
        exportToPPTX('close');
    end
    
end



