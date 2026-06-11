document.addEventListener('DOMContentLoaded', () => {
    const form = document.getElementById('design-form');
    const submitBtn = document.getElementById('submit-btn');
    const btnText = document.getElementById('btn-text');
    const btnSpinner = document.getElementById('btn-spinner');
    
    const welcomeMessage = document.getElementById('welcome-message');
    const errorMessage = document.getElementById('error-message');
    const resultsContent = document.getElementById('results-content');
    
    const successCount = document.getElementById('success-count');
    const failedCount = document.getElementById('failed-count');
    const successRate = document.getElementById('success-rate');
    
    const resultsTableBody = document.querySelector('#results-table tbody');
    const gelContainer = document.getElementById('gel-container');
    const failedSnpsContainer = document.getElementById('failed-snps-container');
    const failedSnpsList = document.getElementById('failed-snps-list');
    
    const downloadBtn = document.getElementById('download-btn');
    
    let currentResults = [];

    form.addEventListener('submit', async (e) => {
        e.preventDefault();
        
        // Reset UI
        errorMessage.style.display = 'none';
        welcomeMessage.style.display = 'none';
        resultsContent.style.display = 'none';
        submitBtn.disabled = true;
        btnText.style.display = 'none';
        btnSpinner.style.display = 'block';
        
        const formData = new FormData(form);
        
        try {
            const response = await fetch('/api/design', {
                method: 'POST',
                body: formData
            });
            
            const data = await response.json();
            
            if (!response.ok) {
                throw new Error(data.detail || 'An error occurred during processing.');
            }
            
            displayResults(data);
            
        } catch (error) {
            errorMessage.textContent = error.message;
            errorMessage.style.display = 'block';
            welcomeMessage.style.display = 'block';
        } finally {
            submitBtn.disabled = false;
            btnText.style.display = 'block';
            btnSpinner.style.display = 'none';
        }
    });

    function displayResults(data) {
        currentResults = data.results;
        
        // Update summary cards
        successCount.textContent = data.success_count;
        failedCount.textContent = data.failed_snps.length;
        
        const rate = data.total_count > 0 ? ((data.success_count / data.total_count) * 100).toFixed(1) : 0;
        successRate.textContent = `${rate}%`;
        
        // Populate table
        resultsTableBody.innerHTML = '';
        data.results.forEach(res => {
            const tr = document.createElement('tr');
            tr.innerHTML = `
                <td>${res['Seq ID']}</td>
                <td>${res['SNP Position']}</td>
                <td>${res['Ref/Alt']}</td>
                <td>${res['Target Tm']}</td>
                <td style="font-family: monospace; font-size: 0.85em;">${res['IF Primer']}</td>
                <td>${res['IF Tm']}</td>
                <td style="font-family: monospace; font-size: 0.85em;">${res['IR Primer']}</td>
                <td>${res['IR Tm']}</td>
            `;
            resultsTableBody.appendChild(tr);
        });
        
        // Failed SNPs
        if (data.failed_snps.length > 0) {
            failedSnpsList.innerHTML = data.failed_snps.map(err => `<li>${err}</li>`).join('');
            failedSnpsContainer.style.display = 'block';
        } else {
            failedSnpsContainer.style.display = 'none';
        }
        
        // Gel visualization
        if (data.gel_data && data.gel_data.length > 0) {
            drawGel(data.gel_data);
            gelContainer.style.display = 'block';
        } else {
            gelContainer.style.display = 'none';
        }
        
        resultsContent.style.display = 'block';
    }
    
    function drawGel(gelData) {
        const traces = [];
        const numLanes = gelData.length;
        
        // Background shape (grey gel block)
        const shapes = [
            {
                type: 'rect',
                x0: -1, y0: 0, x1: numLanes, y1: 100,
                line: { color: 'black', width: 2 },
                fillcolor: 'rgba(245,245,245,0.95)'
            }
        ];
        
        // Wells
        for (let i = 0; i < numLanes; i++) {
            shapes.push({
                type: 'circle',
                x0: i - 0.2, y0: -5, x1: i + 0.2, y1: -1,
                line: { color: 'black' },
                fillcolor: 'black'
            });
        }
        
        // Get min and max sizes to scale
        let allSizes = [];
        gelData.forEach(d => {
            Object.values(d.amplicon_sizes).forEach(size => {
                if (size !== null) allSizes.push(size);
            });
        });
        
        const maxSize = Math.max(...allSizes, 100);
        const minSize = Math.min(...allSizes, 50);
        const rangeSpan = maxSize - minSize || 1;
        
        const colors = ['#1e40af', '#6b21a8', '#d97706', '#b91c1c'];
        const annotations = [];
        
        gelData.forEach((d, laneIdx) => {
            // Lane label
            annotations.push({
                x: laneIdx, y: 105,
                text: `Lane ${laneIdx+1}<br>${d.seq_id}`,
                showarrow: false,
                font: { size: 10 }
            });
            
            const sizes = Object.values(d.amplicon_sizes).filter(s => s !== null);
            sizes.forEach((size, i) => {
                const position = 100 - ((size - minSize) / rangeSpan) * 90;
                const color = colors[i % colors.length];
                
                traces.push({
                    x: [laneIdx-0.3, laneIdx+0.3, laneIdx+0.3, laneIdx-0.3, laneIdx-0.3],
                    y: [position, position, position+4, position+4, position],
                    fill: 'toself',
                    fillcolor: color,
                    line: { color: 'black', width: 1 },
                    hoverinfo: 'text',
                    text: `Size: ${size} bp`,
                    showlegend: false,
                    type: 'scatter',
                    mode: 'lines'
                });
                
                annotations.push({
                    x: laneIdx, y: position + 2,
                    text: String(size),
                    showarrow: false,
                    font: { color: 'white', size: 9, weight: 'bold' }
                });
            });
        });
        
        const layout = {
            title: { text: "Expected Gel Electrophoresis Results", font: { size: 16 } },
            xaxis: { showticklabels: false, showgrid: false, zeroline: false, range: [-1, numLanes] },
            yaxis: { 
                title: "Fragment Size (bp)", 
                range: [-10, 110], 
                autorange: "reversed",
                gridcolor: "lightgray"
            },
            plot_bgcolor: "white",
            paper_bgcolor: "white",
            height: 500,
            margin: { t: 60, b: 20, l: 60, r: 20 },
            shapes: shapes,
            annotations: annotations
        };
        
        Plotly.newPlot('gel-plot', traces, layout, {responsive: true});
    }
    
    downloadBtn.addEventListener('click', () => {
        if (currentResults.length === 0) return;
        
        const headers = Object.keys(currentResults[0]);
        const csvRows = [];
        
        // Add header
        csvRows.push(headers.join(','));
        
        // Add rows
        currentResults.forEach(res => {
            const values = headers.map(header => {
                const val = res[header];
                // Handle commas in values
                return `"${val}"`;
            });
            csvRows.push(values.join(','));
        });
        
        const csvString = csvRows.join('\n');
        const blob = new Blob([csvString], { type: 'text/csv' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = 'tetradyme_results.csv';
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    });
});
