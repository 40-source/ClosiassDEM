<?php

error_reporting(E_ALL);
ini_set("display_errors", 1);
ini_set('precision', 3);

function quatres_facteur($eff,$TR){
    $deformation = (293/$TR)**(1/2); //effet thermique sur la section efficace
    $epsilon=$_GET['epsilon']; // facteur correctif de fissions rapides
    $p=$_GET['p']; //facteur antitrappe
    $v5=$_GET['v5']*$deformation;
    $a5=$_GET['a5']*$deformation; // barns
    $a8=$_GET['a8']*$deformation; // barns
    $f5=$_GET['f5']*$deformation; // 582 barns
    $r=$_GET['r']; // rapport modérateur
    $ra=$_GET['ra']; // rapport absorbeur
    $e=$_GET['e']; //%
    $amod=$_GET['amod']; // barns
    $abarre=$_GET['abarre']; // barns
    $fn=($v5*$f5*$e)/($e*$a5+(1-$e)*$a8+$r*1*$amod+$ra*1*$abarre*$eff);
    $k_inf=$p*$epsilon*$fn;
    return $k_inf;
}

function PNF(){
    $L2=$_GET['L2']; // cm²
    $h_coeur=$_GET['h_coeur']; // cm
    $rayon_coeur = $_GET['rayon_coeur']; // cm
    $Bg=(3.14/$h_coeur)**2+(2.405/$rayon_coeur)**2;
    $PNF=1/(1+$L2*$Bg**2);
    return $PNF;
}

function efficient($k_inf,$PNF,$TR){
    $doppler = ($_GET['doppler']+45/1000*$_GET['Bore']+15)*10**-5; // -pcm/K
    $k_eff=$k_inf*$PNF;
    $k_eff = 1/(1/$k_eff - $doppler*($TR-293)); //effet thermique Doppler
    return $k_eff;
}

function inverse($PNF,$TR){
    $deformation = (293/$TR)**(1/2); //effet thermique sur la section efficace
    $epsilon=$_GET['epsilon']; // facteur correctif de fissions rapides
    $doppler = ($_GET['doppler']+45/1000*$_GET['Bore']+15)*10**-5; // -pcm/K
    $p=$_GET['p']; //facteur antitrappe
    $v5=$_GET['v5']*$deformation;
    $a5=$_GET['a5']*$deformation; // barns
    $a8=$_GET['a8']*$deformation; // barns
    $f5=$_GET['f5']*$deformation; // 582 barns
    $r=$_GET['r']; // rapport modérateur
    $ra=$_GET['ra']; // rapport absorbeur
    $e=$_GET['e']; //%
    $amod=$_GET['amod']; // barns
    $abarre=$_GET['abarre']; // barns

    $vrai_1 = 1/(1 + $doppler*($TR-293));
    $iv=$vrai_1/$PNF;
    $iv=$iv/($p*$epsilon);
    $eff=(($v5*$f5*$e)/$iv-($e*$a5+(1-$e)*$a8+$r*1*$amod))/($ra*1*$abarre);
    return $eff;
}

function reactivité($k_eff){
    $p=($k_eff-1)/$k_eff;
    return $p;
}

function phi($phi,$octave){
    $t=$_GET['t'];
    $phi=$phi*exp($octave*$t);
    return $phi;
}

function puissance_volumique($phi){
    $somme_f=$_GET['somme_f']; // cm⁻¹ section efficace macroscopique
    $énergie_de_fission=3.3*10**(-11);// Energie dégagée par une fission en J
    return $phi*$somme_f*$énergie_de_fission;
}

function JQ($puissance_volumique) {
    $rayon = $_GET['rayon_coeur']*10**-2; //m
    return $rayon*$puissance_volumique/2;
}

function Tgv($Pcomd,$Tgv,$puissance_volumique){
    $h=$_GET['h_coeur']*10**-2; // m
    $rayon = $_GET['rayon_coeur']*10**-2; //m
    $volume = 3.14*$h*$rayon**2; // m3 volume du cylindre
    $dt = $_GET['t']; // s
    $cp_vap = $_GET['Cp_vap']; // j/kg/K
    $m_secondaire = $_GET['masse_secondaire']; //kg

    $Tgv =  $Tgv + $dt * ($volume * $puissance_volumique - $Pcomd)/($m_secondaire * $cp_vap);
    return $Tgv;
}

function puissanceGV($TR,$Tgv,$puissance_volumique){
    $h=$_GET['h_coeur']*10**-2; // m
    $rayon = $_GET['rayon_coeur']*10**-2; //m
    $volume = 3.14*$h*$rayon**2; // m3 volume du cylindre
    $thermostat = $_GET['thermostat']; // K (20°C)
    $dt = $_GET['t']; // s
    $cp_vap = $_GET['Cp_vap']; // j/kg/K
    $m_secondaire = $_GET['masse_secondaire']; //kg
    $Pcomd_lim = $_GET['Plim']*10**6; // MW

    // Automate
    $dT_GV = $_GET['gradiant_GV']*($TR-$thermostat)/(591-$thermostat); //K
    $Tgv_consigne = $TR - $dT_GV; //K
    $Pcomd = $volume * $puissance_volumique - ($Tgv_consigne - $Tgv)/$dt * $m_secondaire * $cp_vap;
    //fin Automate

    if($Pcomd<0){ // hypothèse : Le cricuit de refroidissement ne peut pas chauffer le GV
        return 0;
    }elseif($Pcomd > $Pcomd_lim){ // Le circuit de refroidissement a une puissance maximale qui, même si l'automate le désire, ne peut être dépassée
        return $Pcomd_lim;
    }else{
        return $Pcomd;
    }
}

function TR($TR,$Tgv,$puissance_volumique) {
    $h=$_GET['h_coeur']*10**-2; // m
    $rayon = $_GET['rayon_coeur']*10**-2; //m
    $volume = 3.14*$h*$rayon**2; // m3 volume du cylindre
    $thermostat = $_GET['thermostat']; // K (20°C)
    $cp = $_GET['Cp']; //5829 j/kg/K
    $Débit_massique = $_GET['Dm']; // kg/s
    $masse = $_GET['masse_primaire']; // kg
    $dt = $_GET['t']; // s

    $lambda =  $_GET['lambda'];
    $surface_gv = $_GET['surface_gv']; // m2
    $e = $_GET['ep']; // m
    $Débit_massique_gv = $_GET['Dm_gv']; //kg/s

    $Tf = ($TR - $Tgv)*exp(-($lambda*$surface_gv)/($cp*$Débit_massique_gv*$e)) + $Tgv;
    $dT_GV = $TR - $Tf;

    $thermo_forcage = ($puissance_volumique * $volume)/($cp * $masse) - $Débit_massique * $dT_GV / $masse;
    return $TR + $thermo_forcage*$dt;
}

function puissanceTh($JQ,$surface_combustible){
    return $JQ*$surface_combustible;
}

function octave($k_eff){
    $tau=$_GET['tau']; //s
    $beta=$_GET['beta']; //sans dimension
    $l=$_GET['l']; //s
    $l_prec=$l+$beta*$tau;
    return ($k_eff-1)/$l_prec;
}

function puissanceTurbine($puissanceGV){
    $débit = $_GET['débit_massique'];
    $h2=$_GET['h2']; //kj/kg
    $h4=$_GET['h4']; //kj/kg
    $h3=$puissanceGV/$débit*10**-3+$h2; //kj/kg
    $wt=$h4-$h3; //kj/kg
    if($wt*$débit*10**3<0){
        return $wt*$débit*10**3; //W
    }else{
        return 0;
    }
}

function RED($value) {
    $minValue = 293;
    $maxValue = 800;
    $proportion = ($value - $minValue) / ($maxValue - $minValue);
    $proportion = 1 - $proportion;
    $red = round($proportion * (255 - 255/6)); // orange à rouge
    $green = round($proportion * (255/2) + 255/2); // rouge à jaune
    $blue = round($proportion * (255/3) + 2 * 255/3); // jaune à pourpre
    return sprintf("#%02x%02x%02x", $blue, $red, $red);
}

function Yellow($value) {
    $minValue = 0;
    $maxValue = 130;
    $proportion = ($value - $minValue) / ($maxValue - $minValue);
    $proportion = 1 - $proportion;
    $red = round($proportion * (255 - 255/6)); // orange à rouge
    $green = round($proportion * (255/2) + 255/2); // rouge à jaune
    $blue = round($proportion * (255/3) + 2 * 255/3); // jaune à pourpre
    return sprintf("#%02x%02x%02x",$green , $blue, $red);
}

function Yellowa($value) {
    $minValue = 0;
    $maxValue = 5000;
    $proportion = ($value - $minValue) / ($maxValue - $minValue);
    $proportion = 1 - $proportion;
    $red = round($proportion * (255 - 255/6)); // orange à rouge
    $green = round($proportion * (255/2) + 255/2); // rouge à jaune
    $blue = round($proportion * (255/3) + 2 * 255/3); // jaune à pourpre
    return sprintf("#%02x%02x%02x",$green , $blue, $red);
}

function Yellowb($value) {
    $minValue = 0;
    $maxValue = 2900;
    if($value>0){
        $proportion = ($value - $minValue) / ($maxValue - $minValue);
        $proportion = 1 - $proportion;
        $red = round($proportion * (255 - 255/6)); // orange à rouge
        $green = round($proportion * (255/2) + 255/2); // rouge à jaune
        $blue = round($proportion * (255/3) + 2 * 255/3); // jaune à pourpre
        return sprintf("#%02x%02x%02x",$green , $blue, $red);
    }else{
        $proportion = (-$value - $minValue) / ($maxValue - $minValue);
        $proportion = 1 - $proportion;
        $red = round($proportion * (255 - 255/6)); // orange à rouge
        $green = round($proportion * (255/2) + 255/2); // rouge à jaune
        $blue = round($proportion * (255/3) + 2 * 255/3); // jaune à pourpre
        return sprintf("#%02x%02x%02x",$red , $blue, $green);
    }
    
}

function Yellowc($value) {
    $minValue = 0;
    $maxValue = 100;
    $proportion = ($value - $minValue) / ($maxValue - $minValue);
    $proportion = 1 - $proportion;
    $red = round($proportion * (255 - 255/6)); // orange à rouge
    $green = round($proportion * (255/2) + 255/2); // rouge à jaune
    $blue = round($proportion * (255/3) + 2 * 255/3); // jaune à pourpre
    return sprintf("#%02x%02x%02x",$green , $blue, $red);
}

function Purple($value) {
    $minValue = 10**11;
    $maxValue = 10**15;
    $proportion = ($value - $minValue) / ($maxValue - $minValue);
    $proportion = 1 - $proportion;
    $red = round($proportion * (255 - 255/6)); // orange à rouge
    $green = round($proportion * (255/2) + 255/2); // rouge à jaune
    $blue = round($proportion * (255/3) + 2 * 255/3); // jaune à pourpre
    return sprintf("#%02x%02x%02x",$red , $green, $blue);
}

function Purplea($value) {
    $minValue = 0;
    $maxValue = 1;
    $proportion = ($value - $minValue) / ($maxValue - $minValue);
    $proportion = 1 - $proportion;
    $red = round($proportion * (255 - 255/6)); // orange à rouge
    $green = round($proportion * (255/2) + 255/2); // rouge à jaune
    $blue = round($proportion * (255/3) + 2 * 255/3); // jaune à pourpre
    return sprintf("#%02x%02x%02x",$red , $green, $blue);
}

function Green($value) {
    $red = max(0, min(255, round(255 * abs($value - 1))));
    $green = max(0, min(255, round(255 * (1 - abs($value - 1)))));
    $blue = 0;
    return sprintf("#%02x%02x%02x", $red, $green, $blue);
}

$i=0;
$t=$_GET['t'];
$total=$_GET['total'];
$phi=$_GET['phi'];
$PNF=PNF();
$surface_combustible=$_GET['h_coeur']*10**-2*($_GET['rayon_coeur']*10**-2)*2*3.14;
$fuite=2; //2% de fuite

$amplitude=$_GET['amplitude'];
$T_consigne=$_GET['T_consigne'];
$pv_consigne  =$_GET['pv_consigne'];
$sécurité=$_GET['sécurité'];

$TR=293; // K
$Tgv=293; // K
$Tsat = 561.15; // K

$eff_ad = 0;
$eff_pu = 0;
$amod=$_GET['amod']*(1-$_GET['Bore']/10**6)+$_GET['ab']*$_GET['Bore']/10**6;
$_GET['amod']=$amod;
$eff=inverse($PNF,$TR);

echo  "(Borication tel que a<sub>mod</sub> = ". number_format((float)$amod, 2, '.', '')." cm⁻¹)";

echo "<table>";

error_reporting(21495);
ini_set("display_errors", 0);

echo "<tr><td>";
echo "<td>Temps<td>Enfoncement des barres<td>Coefficient de reproduction infini<td>Coefficient de reproduction efficace<td>Pulsation principale<td>Octavemètre<td>Réactivité<td>Flux surfacique de neutron<td>Puissance volumique<td>Température primaire<td>Température secondaire<td>Puissance réacteur<td>Puissance secondaire (gv/comd)<td>Puissance à la turbine<td>Rendement<td>Consigne de l'automate";
$SCRAM=0;

while($i<$total){
    $k_inf=quatres_facteur($eff,$TR);
    $k_eff=efficient($k_inf,$PNF,$TR);
    $octave=octave($k_eff);
    $p=reactivité($k_eff);
    $phi=phi($phi,$octave);
    $pv=puissance_volumique($phi);
    $JQ=JQ($pv*10**6);
    $TR=TR($TR,$Tgv,$pv*10**6);
    $Pgv=puissanceGV($TR,$Tgv,$pv*10**6);
    $Tgv=Tgv($Pgv,$Tgv,$pv*10**6);
    $Pth=puissanceTh($JQ,$surface_combustible);
    $Pt=puissanceTurbine($Pgv);
    $r=intval(-1*$Pt/$Pgv*100)-$fuite;

    if($r<0){
        $r=0;
    }
    if($r>100){
        $r="∅";
    }

    //Asservissement
    if($pv<$pv_consigne and $TR<$T_consigne){
        if($sécurité>0){
            if($eff-$amplitude*($pv_consigne-$pv)/$pv>0){
                $eff_ad=$eff_ad-$amplitude*($pv_consigne-$pv)/$pv;
                $eff = inverse($PNF,$TR) + $eff_ad;
                $eff_pu = 0;
            }else{
                $eff=0;
            }
        }else{
            $eff=inverse($PNF,$TR)-$amplitude*($pv_consigne-$pv)/$pv;
        }
    }elseif ($k_eff>0 or $phi>10**15) {
        if($sécurité>0){
            if($eff+$amplitude*($pv-$pv_consigne)/$pv<100){
                $eff_pu=$eff_pu+$amplitude*($pv-$pv_consigne)/$pv;
                $eff = inverse($PNF,$TR) + $eff_pu;
                $eff_ad = 0;
            }else{
                $eff=100;
            }
        }else{
            $eff=inverse($PNF,$TR)+$amplitude*($TR-$T_consigne)/$TR;
        }
    }

    $date[$i][0]="ON";
    //Consigne
    if($Pt>=0){
        $consigne="Démarrage ☢️";
        $indicteur_turbine=0;
        $color_consigne="Yellow";
        $date[$i][0]="OFF";
        $date[$i][6]="DEMN";
    }elseif($TR<$T_consigne-10){
        $consigne="Montée en puissance ⚠️";
        $color_consigne="Red";
        $date[$i][6]="MOPU";
    }elseif($TR>$T_consigne+10){
        $consigne="Baisse de puissance ⚠️";
        $color_consigne="Red";
        $date[$i][6]="BAPU";
    }else{
        $consigne="STABLE 👍";
        $color_consigne="Green";
        $date[$i][6]="STAB";
    }


    //Style :
    if($indicteur_turbine==1){
        $consigne="Turbine connectée ⚡️";
        $color_consigne="Purple";
        $indicteur_turbine=2;
    }elseif($indicteur_turbine!=2){
        $indicteur_turbine=1;
    }
    
    if($p>0){
        $color_réactivité="red";
    }else{
        $color_réactivité="green";
    }

    $date[$i][1]=round($k_eff,3);
    $date[$i][2]=intval($TR);
    $date[$i][3]=intval(-$Pt*10**-6)." MW";
    
    if($sécurité<2){
        if($TR>900 or $k_eff>1.7 or -$Pt*10**-6>2000 or $phi>10**16 or $SCRAM==1){
            $consigne="SCRAM";
            $eff=inverse($PNF,$TR)*2;
            $date[$i][3]="SCRAM ⚠";
            $date[$i][6]="SCRAM";
            $SCRAM=1;
            $color_consigne="Red";
        }
    }
    
    $date[$i][4]=$eff*100;
    $date[$i][5]=round(log10($phi),2);
    $date[$i][7]=round($octave*60/0.69,3);
    $date[$i][8]=round($p*10**5,3);
    $date[$i][9]=round(log10($phi),2);
    $date[$i][10]=intval($Pth*10**-6)." MW";
    $date[$i][11]=intval($Tgv);

    $color_T=RED(intval($TR));
    $color_T1=RED(intval($Tgv));
    $color_Pv=Yellow(intval($pv));
    $color_Phi=Purple(intval($phi));
    $color_keff=Green($k_eff);
    $color_inf=Green($k_inf);

    $color_Pgv=Yellowa(intval($Pgv*10**-6));
    $color_Pt=Yellowb(intval(-$Pt*10**-6));
    $color_r=Yellowc($r);

    $color_amod=Purplea($eff);

    echo "<tr><td>";
    echo "<td> t = ". ($t*$i)." s";
    echo "<td style='background-color: ".$color_amod.";'> enf<sub>barres</sub> = ". round($eff*100,2)." %";
    echo "<td style='background-color: ".$color_inf.";'> k<sub>∞</sub> = ". round($k_inf,3);
    echo "<td style='background-color: ".$color_keff.";'> k<sub>eff</sub> = ". round($k_eff,3);
    echo "<td style='background-color: ".$color_keff.";'> ω = ". round($octave,3) ." s⁻¹";
    echo "<td style='background-color: ".$color_keff.";'> Ω = ". round($octave*60/0.69,3). " octave.mn⁻¹";
    echo "<td style='background-color: ".$color_réactivité.";'> ρ = ". round($p,5);
    echo "<td style='background-color: ".$color_Phi.";'> φ = ".intval($phi)." n cm⁻² s⁻¹";
    echo "<br>log<sub>10</sub> = ".(round(log10($phi),2));
    echo "<td style='background-color: ".$color_Pv.";'> 𝛿 = ".intval($pv)." MW/m³";
    echo "<td style='background-color: ".$color_T.";'> T<sub>pr</sub> = ".intval($TR)." K";
    if($Tgv<$Tsat){
        echo "<td style='background-color: ".$color_T1.";'> T<sub>gv</sub> = ".intval($Tgv)." K";
    }else{
        echo "<td style='background-color: red;'> T<sub>gv</sub> = ".intval($Tgv)." K ♨️ Saturation";
    }
    echo "<td style='background-color: ".$color_Pv.";'> P<sub>th</sub> = ".intval($Pth*10**-6)." MW";
    echo "<td style='background-color: ".$color_Pgv.";'> P<sub>gv</sub> = ".intval($Pgv*10**-6)." MW";
    echo "<td style='background-color: ".$color_Pt.";'> P<sub>t</sub> = ".intval(-$Pt*10**-6)." MW";
    echo "<td style='background-color: ".$color_r.";'> η = ".$r." %";
    echo "<td style='background-color: ".$color_consigne.";'> ".$consigne;
    echo "<td><tr>";
    

    $i=$i+1;
}
$id=0;

while($TR>294 and $i<$total*2){
    $k_inf=quatres_facteur($eff,$TR);
    $k_eff=efficient($k_inf,$PNF,$TR);
    $octave=octave($k_eff);
    $p=reactivité($k_eff);
    $phi=phi($phi,$octave);
    $pv=puissance_volumique($phi);
    $JQ=JQ($pv*10**6);
    $TR=TR($TR,$Tgv,$pv*10**6);
    $Pgv=puissanceGV($TR,$Tgv,$pv*10**6);
    $Tgv=Tgv($Pgv,$Tgv,$pv*10**6);
    $Pth=puissanceTh($JQ,$surface_combustible);
    $Pt=0;
    $r=intval(-1*$Pt/$Pgv*100)-$fuite;

    if($r<0){
        $r=0;
    }
    if($r>100){
        $r="∅";
    }

    //Asservissement
    if($pv>2){
        $eff=inverse($PNF,$TR)+$amplitude*($pv-2)/$pv;
    }

    $date[$i][0]="ON";
    $date[$i][6]="BAPU";
    //Consigne
    if($Pt>=0){
        $consigne="Arret ☢️";
        $indicteur_turbine=0;
        $color_consigne="Yellow";
        $date[$i][0]="OFF";
        $date[$i][6]="SHUT";
    }elseif($TR<$T_consigne-10){
        $consigne="Baisse de puissance ⚠️";
        $color_consigne="Red";
        $date[$i][6]="BAPU";
    }

    //Style :
    if($indicteur_turbine==0 and $id==0){
        $consigne="Turbine déconnectée ⚡️";
        $color_consigne="Purple";
        $id=1;
        $date[$i][6]="BAPU";
    }

    //Style :
    if($indicteur_turbine==1){
        $consigne="Turbine connectée ⚡️";
        $color_consigne="Purple";
        $indicteur_turbine=2;
    }elseif($indicteur_turbine!=2){
        $indicteur_turbine=1;
    }
    
    if($p>0){
        $color_réactivité="red";
    }else{
        $color_réactivité="green";
    }

    $color_T=RED(intval($TR));
    $color_T1=RED(intval($Tgv));
    $color_Pv=Yellow(intval($pv));
    $color_Phi=Purple(intval($phi));
    $color_keff=Green($k_eff);
    $color_inf=Green($k_inf);

    $color_Pgv=Yellowa(intval($Pgv*10**-6));
    $color_Pt=Yellowb(intval(-$Pt*10**-6));
    $color_r=Yellowc($r);

    $color_eff=Purplea($eff);

    echo "<tr><td>";
    echo "<td> t = ". ($t*$i)." s";
    echo "<td style='background-color: ".$color_amod.";'> enf<sub>barres</sub> = ". round($eff*100,2)." %";
    echo "<td style='background-color: ".$color_inf.";'> k<sub>∞</sub> = ". round($k_inf,3);
    echo "<td style='background-color: ".$color_keff.";'> k<sub>eff</sub> = ". round($k_eff,3);
    echo "<td style='background-color: ".$color_keff.";'> ω = ". round($octave,3) ." s⁻¹";
    echo "<td style='background-color: ".$color_keff.";'> Ω = ". round($octave*60/0.69,3). " octave.mn⁻¹";
    echo "<td style='background-color: ".$color_réactivité.";'> ρ = ". round($p,5);
    echo "<td style='background-color: ".$color_Phi.";'> φ = ".intval($phi)." n cm⁻² s⁻¹";
    echo "<br>log<sub>10</sub> = ".(round(log10($phi),2));
    echo "<td style='background-color: ".$color_Pv.";'> 𝛿 = ".intval($pv)." MW/cm³";
    echo "<td style='background-color: ".$color_T.";'> T<sub>pr</sub> = ".intval($TR)." K";
    echo "<td style='background-color: ".$color_T1.";'> T<sub>gv</sub> = ".intval($Tgv)." K";
    echo "<td style='background-color: ".$color_Pv.";'> P<sub>th</sub> = ".intval($Pth*10**-6)." MW";
    echo "<td style='background-color: ".$color_Pgv.";'> P<sub>comd</sub> = ".intval($Pgv*10**-6)." MW";
    echo "<td style='background-color: ".$color_Pt.";'> P<sub>t</sub> = ".intval(-$Pt*10**-6)." MW";
    echo "<td style='background-color: ".$color_r.";'> η = ".$r." %";
    echo "<td style='background-color: ".$color_consigne.";'> ".$consigne;
    echo "<td><tr>";

    $date[$i][1]=round($k_eff,3);
    $date[$i][2]=intval($TR);
    $date[$i][3]=intval(-$Pt*10**-6)." MW";
    $date[$i][4]=$eff*100;
    $date[$i][5]=round(log10($phi),2);
    $date[$i][7]=round($octave*60/0.69,3);
    $date[$i][8]=round($p*10**5,3);
    $date[$i][9]=round(log10($phi),2);
    $date[$i][10]=intval($Pth*10**-6)." MW";
    $date[$i][11]=intval($Tgv);

    $i=$i+1;
}
echo "<tr><td>";
echo "<td> END</table>";
$j=json_encode($date);
?>
<script type="text/javascript">
    var cpt = 0;
    var url = "img5/aff.php";
    var lim=<?php echo $i;?>;

    function ajax()
    {
        j=<?php echo $j;?>;
        $.ajax({url: url, success: function(result)
            {
                $("#chat").empty();
                $("#chat").html(result);
                $("#chat").ready(function() {
                    if (cpt<lim){
                        cpt=cpt+1;
                        url = "img5/aff.php?i="+cpt+"&T="+j[cpt][2]+"&K="+j[cpt][1]+"&P="+j[cpt][3]+"&mode="+j[cpt][0]+"&eff="+j[cpt][4]+"&phi="+j[cpt][5]+"&consigne="+j[cpt][6]+"&octave="+j[cpt][7]+"&réactivité="+j[cpt][8]+"&flux="+j[cpt][9]+"&Pth="+j[cpt][10]+"&Tgv="+j[cpt][11];
                        setTimeout(ajax, 100);
                    }
                });
            }});
    }

    ajax();
</script>