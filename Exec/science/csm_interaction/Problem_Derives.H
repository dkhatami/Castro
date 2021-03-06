derive_lst.add("cd_mask",IndexType::TheCellType(),1,ca_dercd_mask,the_same_box);
derive_lst.addComponent("cd_mask",desc_lst,State_Type,Density,1);
derive_lst.addComponent("cd_mask",desc_lst,State_Type,Xmom,3);
derive_lst.addComponent("cd_mask",desc_lst,State_Type,Eint,1);
derive_lst.addComponent("cd_mask",desc_lst,State_Type,Shock,1);
derive_lst.addComponent("cd_mask",desc_lst,State_Type,FirstAux,2);

derive_lst.add("csm_mask",IndexType::TheCellType(),1,ca_dercsm_mask,the_same_box);
derive_lst.addComponent("csm_mask",desc_lst,State_Type,Density,1);
derive_lst.addComponent("csm_mask",desc_lst,State_Type,FirstAux+1,1);
derive_lst.addComponent("csm_mask",desc_lst,State_Type,FirstAux+2,1);

derive_lst.add("dtau",IndexType::TheCellType(),1,ca_derdtau,the_same_box);
derive_lst.addComponent("dtau",desc_lst,State_Type,Density,1);

derive_lst.add("tau_r",IndexType::TheCellType(),1,ca_dertau_r,the_same_box);
