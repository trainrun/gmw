head = """<!DOCTYPE html>
<html>
<head>
    <title>GMW Debruijn Graph</title>
    <script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
    <link href="https://unpkg.com/vis-network/styles/vis-network.css" rel="stylesheet" type="text/css" />
</head>
<body>
    <h1 style="text-align: center;"></>GMW Debruijn Graph</h1>
    <h4 style="text-align: center;"></>Each node in the graph represent a contig, the length of the contig are signed below the nodes. Pink for target contigs, Blue for contaminate congtis, and Grey for unknown contigs.</h4>
    <div id="gmw" style="width: 100%; height: 800px; margin: 0 auto;"></div>

    <script>
        var jsonData = 
"""

middle = """;
        var options = 
"""
tail = """;
        var nodes = new vis.DataSet(jsonData.nodes);
        var edges = new vis.DataSet(jsonData.edges);

        var data = {
            nodes: nodes,
            edges: edges
        };
        var network = new vis.Network(document.getElementById('gmw'), data, options);
    </script>
</body>
</html>
"""