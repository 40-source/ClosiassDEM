function fichierExiste(url, callback) {
    var xhr = new XMLHttpRequest();
    xhr.onreadystatechange = function() {
        if (xhr.readyState == 4) {
            callback(xhr.status == 200);
        }
    };
    xhr.open('HEAD', url, true);
    xhr.send();
}

function range(start, end) {
    var ans = [];
    for (let i = start; i <= end; i++) {
        ans.push(i);
    }
    return ans;
}

let combx = []
let comby = []
let combz = []

d3.csv("./"+nom+ "/comb.csv", function(err, rows) {
    function unpack(rows, key) {
        return rows.map(function(row) { return row[key]; });
    }
    combx = unpack(rows, 'x')
    comby = unpack(rows, 'y')
    combz = unpack(rows, 'z')

});

function affichage(url){
d3.csv(url, function(err, rows){

    function unpack(rows, key) {
        return rows.map(function(row) { return row[key]; });
    }

    var data = [{
        x: unpack(rows, 'x'),
        y: unpack(rows, 'y'),
        z: unpack(rows, 'z'),
        name: "Neutrons",
        mode: 'markers',
        type: 'scatter3d',
        marker: {
          color: 'rgb(23, 190, 207)',
          size: 1
        }
    }];

    var data0 = [{
        x: combx,
        y: comby,
        z: combz,
        name: "Comb",
        mode: 'markers',
        type: 'scatter3d',
        marker: {
          color: 'rgb(58, 255, 0)',
          size: 1
        }
    }];

    ne.push(rows.length)
    if (!isNaN(parseFloat(rows.length / ne[ne.length - 2]))) {
        total = total + rows.length / ne[ne.length - 2];
    }
    k.push(total/i)

    var trace = {
        x: range(0,i),
        y: ne,
        type: 'scatter'
    };
    data1 = [trace]

    var titleN = {
        height: 580,
        title:'Variation de neutrons'
    };

    Plotly.newPlot('n', data1, titleN);

    var trace1 = {
        x: range(0,i),
        y: k,
        type: 'scatter'
    };
    data2 = [trace1]

    var titleK = {
        height: 580,
        title:'k_eff'
    };
    Plotly.newPlot('k', data2, titleK);

    var layout = {
        autosize: true,
        height: 580,
        scene: {
            aspectratio: {
                x: 1,
                y: 1,
                z: 1
            },
            camera: {
                center: {
                    x: 0,
                    y: 0,
                    z: 0
                },
                eye: {
                    x: 0,
                    y: 0,
                    z: 2
                },
                up: {
                    x: 0,
                    y: 0,
                    z: 0
                }
            },
            xaxis: {
                type: 'linear',
                zeroline: false,
                range: [-30, 30]
            },
            yaxis: {
                type: 'linear',
                zeroline: false,
                range: [-30, 30]
            },
            zaxis: {
                type: 'linear',
                zeroline: false,
                range: [-30, 30]
            }
        },
        title: 'Closiass DEM'
    };

    Plotly.newPlot('aff', data, layout);
    Plotly.addTraces('aff', data0);

});
}

let i = 1;
let ne = [];
let k = [];
let total = 1;
let urlDuFichier = "./" + nom + "/" + i + ".csv";
let urlEnTexte = urlDuFichier.toString();
let isRunning = true;

function verifierEtAfficher() {
    if (!isRunning) {
        setTimeout(verifierEtAfficher, 2);
        return;
    }
    fichierExiste(urlEnTexte, function(existe) {
        if (existe) {
            i++;
            affichage(urlEnTexte)
            urlDuFichier = "./" + nom + "/" + i + ".csv";
            urlEnTexte = urlDuFichier.toString();
        }
        setTimeout(verifierEtAfficher, 20);
    });
}

function toggleExecution() {
    isRunning = !isRunning;
}

verifierEtAfficher()