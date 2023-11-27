<!DOCTYPE html>
<html lang="fr">
<link rel="icon" href="../logo.png" type="png">
<head>
  <meta charset="UTF-8">
  <title>Movie C</title>
</head>
<body>
<!-- partial:index.partial.html -->
<head>
	<!-- Load plotly.js into the DOM -->
	<script src='plotly-2.27.0.min.js'></script>
	<script src='d3.min.js'></script>
</head>
<center>
	<body>
		<button onclick="toggleExecution()">▶️/⏸</button><button onclick="location.reload()">🔄</button>
		<div id='aff'></div>
		<div id='k'></div>
		<div id='n'></div>
	</body>
</center>
<!-- partial -->
<script>
	var nom = "<?php if(isset($_GET['nom'])){echo $_GET['nom'];}else{echo "default";} ?>"
</script>
<script src="./script.js"></script>

</body>
</html>
