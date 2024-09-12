%Pour chacun des cas faire attention que n>21 pour que la fenetre ne soit pas confoncue avec le mur


%Pour la partie 1 :

  %Pour le test 1 sans chauffage :
    %executer fichiers Chambre1_Simu_Statique_SC et Chambre2_Simu_Statique_SC
      %Pour la chambre 1 executez par exemmple la commande suivante :
        Chambre1SimuStatiqueSC(20,20,30)
      %Pour la chambre 2 executez la commande suivante :
        Chambre2SimuStatiqueSC(20,20,30)

     

  %Pour le test 2 sans chauffage :
    %executer fichiers Chambre1_Simu_Statique_SC et Chambre2_Simu_Statique_SC
        %Pour la chambre 1 executez par exemmple la commande suivante :
          Chambre1SimuStatiqueSC(-15,10,30)
        %Pour la chambre 2 executez la commande suivante :
          Chambre2SimuStatiqueSC(-15,10,30)



  %Pour le test 3 avec chauffage :
    %executer fichiers Chambre1_Simu_Statique_AC et Chambre2_Simu_Statique_AC
      %Pour la chambre 1 executez par exemmple la commande suivante :
        Chambre1SimuStatiqueAC(-15,10,140,30)
      %Pour la chambre 2 executez la commande suivante :
        Chambre2SimuStatiqueAC(-15,10,500,30)
      
      
      
%Pour la partie 2 :

  %Pour le test 1 (avec chauffage):
    %executer fichiers Chambre1_Simu_Instationnaire et Chambre2_Simu_Instationnaire
      %Pour la chambre 1 executez par exemmple la commande suivante :
        Chambre1SimuInstationnaire(-15,10,140,30)
      %Pour la chambre 2 executez la commande suivante :
        Chambre2SimuInstationnaire(-15,10,500,30)

  
  %Pour le test 2 (avec clim):
    %executer fichiers Chambre1_Simu_Instationnaire et Chambre2_Simu_Instationnaire
      %Pour la chambre 1 executez par exemmple la commande suivante :
        Chambre1SimuInstationnaire(30,25,-200,30)
      %Pour la chambre 2 executez la commande suivante :
        Chambre2SimuInstationnaire(30,25,-500,30)



