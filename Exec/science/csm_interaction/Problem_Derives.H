derive_lst.add("fs_mask",IndexType::TheCellType(),1,ca_derfs_mask,the_same_box);
derive_lst.addComponent("fs_mask",desc_lst,State_Type,URHO,1);
derive_lst.addComponent("fs_mask",desc_lst,State_Type,UFX+3,1);
derive_lst.addComponent("fs_mask",desc_lst,State_Type,Shock,1);

derive_lst.add("csm_mask",IndexType::TheCellType(),1,ca_dercsm_mask,the_same_box);
derive_lst.addComponent("csm_mask",desc_lst,State_Type,URHO,1);
derive_lst.addComponent("csm_mask",desc_lst,State_Type,UFX+3,1);


derive_lst.add("ej_mask",IndexType::TheCellType(),1,ca_derej_mask,the_same_box);
derive_lst.addComponent("ej_mask",desc_lst,State_Type,URHO,1);
derive_lst.addComponent("ej_mask",desc_lst,State_Type,UFX+2,1);

derive_lst.add("dtau",IndexType::TheCellType(),1,ca_derdtau,the_same_box);
derive_lst.addComponent("dtau",desc_lst,State_Type,URHO,1);
derive_lst.addComponent("dtau",desc_lst,State_Type,UTEMP,1);
derive_lst.addComponent("dtau",desc_lst,State_Type,UFX,1);

derive_lst.add("dtau_abs",IndexType::TheCellType(),1,ca_derdtauabs,the_same_box);
derive_lst.addComponent("dtau_abs",desc_lst,State_Type,URHO,1);
derive_lst.addComponent("dtau_abs",desc_lst,State_Type,UTEMP,1);
derive_lst.addComponent("dtau_abs",desc_lst,State_Type,UFX,1);
