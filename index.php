<html>
<link rel="icon" href="../logo.png" type="png">
<style type="text/css">
p
{
   width: 400px;
   border: 1px solid black;
   text-align: justify;
   padding: 12px;
   margin: 50px; /* Marge extérieure de 50px */
}
.description {
            font-size: 18px;
            color: #666;
            margin-bottom: 20px;
}
.source-img {
            margin: 10px;
            vertical-align: middle;
        }

        .source-list {
            list-style: none;
            padding: 0;
            display: flex;
            justify-content: center;
        }

        .source-list-item {
            margin: 10px;
        }

        .source-name {
            font-size: 16px;
            color: #777;
        }
textarea{ resize:none;}
</style>
    <title>DEM</title><meta charset="UTF-8">
<center>
<FONT face="Helvetica">
<?php

if(isset($_POST['n'])&&isset($_POST['t'])&&isset($_POST['r'])&&isset($_POST['uma'])&&isset($_POST['d'])&&isset($_POST['b'])&&isset($_POST['e'])){
        $n=number_format($_POST['n'],1);
        $t=number_format($_POST['t'],1);
        $r=number_format($_POST['r'],1);
        $rep=number_format($_POST['rep'],1);
        $uma=number_format($_POST['uma'],1);
        $d=number_format($_POST['d'],1);
        $b=number_format($_POST['b'],1);
        $e=number_format($_POST['e'],1);
        $commande="nohup ./DEM $n $t $r $rep $uma $d $b $e > /dev/null 2>&1 &";
        exec($commande);
        $nom=$n.$t.$r.$rep.$uma.$d.$b.$e;
        header("Location: play.php?nom=".$nom);
}

?>
    <h1>Simulateur DEM-V°C, Mesolvarde</h1>
    <h2>Diffusion-Émission-Modération</h2>
    <ul class="source-list">
        <li class="source-list-item">
            <img src="../dem/fr.png" alt="Logo de sources françaises" height="50" class="source-img">
            <div class="source-name">PÉRENNOU Kaourant, Commissariat à l'énergie atomique (CEA)<br>HUTIN Gérémy, École normale supérieure de Lyon (ENS)<br>BARRANGER, Noé, Télécom Saint-Étienne (TSE)<br>aprt. Orano Framatome EDF TechnicAtome</div>
        </li>
        <li class="source-list-item">
            <img src="../dem/sn.png" alt="Logo de sources chinoises" height="50" class="source-img">
            <div class="source-name">WANG Tuan, Advanced Semiconductor Engineering (ASE)<br> aprt. CNNC CGNPC Huaneng group</div>
        </li>
    </ul>
<form action="" method="POST">
<p>
<table>
    <thead>
        <tr>
            <th colspan="2">Valeur pour le modèle</th>
        </tr>
    </thead>
    <tbody>
	<tr><td>Nombre de cycle : </td><td><input type="number" name="n" value="20" max="600" min="2"></td></tr>
	<tr><td>Temps par cycle : </td><td><input type="number" name="t" value="0.1" step="0.05" min="0.05" max="0.5"></td></tr>
    <tr><td>Rayon de l'assemblage: </td><td><input type="number" name="r" value="10" min="1"></td></tr>
    <tr><td>Zones de l'assemblage: </td><td><input type="number" name="rep" value="3" min="1"></td></tr>
    <tr><td>Rayon d'action nucléaire : </td><td><input type="number" name="uma" value="1" min="1"></td></tr>
    <tr><td>Hauteur des barres : </td><td><input type="number" name="d"  value="12" min="1"></td></tr>
    <tr><td>Nombre de barres : </td><td><input type="number" name="b"  value="6" min="1"></td></tr>
    <tr><td>Enrichissement : </td><td><input type="number" name="e"  value="2" step="0.1" min="1"></td></tr>
	<tr><td> </td><td><input type="submit" name="submit" value="Transmettre"></td></tr>
        </tbody>
</table>
</p></form>
</center>
