function   [error1,psnr1,error2,psnr2,itr1,itr,error12,psnr12,error22,psnr22,itr12,itr2,error13,psnr13,error23,psnr23,itr13,itr3,error121,psnr121,error221,psnr221,itr121,itr21,error131,psnr131,error231,psnr231,itr131,itr31,error122,psnr122,error222,psnr222,itr122,itr22]=N


  [C1,C,omega]=InitializationTproductirc(100,100,0.1,0.4);
  [error1,psnr1,error2,psnr2,itr1,itr]=M1(C1,C,omega);
  [C1,C,omega]=InitializationTproductirc(100,100,0.1,0.5);
  [error12,psnr12,error22,psnr22,itr12,itr2]=M1(C1,C,omega);
  [C1,C,omega]=InitializationTproductirc(100,100,0.1,0.6);
  [error13,psnr13,error23,psnr23,itr13,itr3]=M1(C1,C,omega);
   [C1,C,omega]=InitializationTproductirc(100,100,0.2,0.3);
  [error121,psnr121,error221,psnr221,itr121,itr21]=M1(C1,C,omega);
   [C1,C,omega]=InitializationTproductirc(100,100,0.2,0.4);
  [error131,psnr131,error231,psnr231,itr131,itr31]=M1(C1,C,omega);
    [C1,C,omega]=InitializationTproductirc(100,100,0.2,0.5);
  [error122,psnr122,error222,psnr222,itr122,itr22]=M1(C1,C,omega);