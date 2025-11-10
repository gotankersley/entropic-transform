window.showResults = function(data) {	
	var resultsDiv = document.getElementById('results');
	var menu = window.menu;
	var template = 
`<table border="1">
	<tr>
		<th>Type</th>
		<th>Sequence</th>
		<th>Rank</th>
		<th>Entropy</th>
	</tr>
	<tr>
		<td>Random</td>
		<td><div class="seq">${data.rndSeq}</div></td>
		<td><div class="rank">${data.rndRank}</div></td>
		<td class="entropy">${data.rndEntropy}</td>
	</tr>`;
	
	if (menu.nearEntropic) {
		template += 
		`<tr>
			<td>Near Entropic</td>
			<td><div class="seq">${data.nearEntShapedSeq}</div></td>
			<td><div class="rank">${data.nearEntRank}</div></td>
			<td class="entropy">${data.nearEntShapedEntropy}</td>
		</tr>`;
	}
	
	if (menu.nearerEntropic) {
		template += 
		`<tr>
			<td>Nearer Entropic</td>
			<td><div class="seq">${data.entShapedSeq}</div></td>
			<td><div class="rank">${data.entRank}</div></td>
			<td class="entropy">${data.entShapedEntropy}</td>
		</tr>`;
	}
	
	if (menu.nearEntropic) {
		template += 
		`<tr>
			<td>Near Entropic  &nbsp;&nbsp;&nbsp;| Base-A</td>
			<td><div class="seq">${data.nearSeqBaseA}</div></td>
			<td><div class="rank">${data.rndRank}</div></td>
			<td class="entropy">${data.nearEntropyBaseA}</td>
		</tr>`;
	}
		
	if (menu.nearerEntropic) {
		template += 
		`<tr>
			<td>Nearer Entropic | Base-A</td>
			<td><div class="seq">${data.entSeqBaseA}</div></td>
			<td><div class="rank">${data.rndRank}</div></td>
			<td class="entropy">${data.entEntropyBaseA}</td>
		</tr>`;	
	}

	if (menu.bwts) {
		template += 
		`<tr>
			<td>BWTS</td>
			<td><div class="seq">${data.bwtsMtfSeq}</div></td>
			<td><div class="rank">${data.bwtsMtfRank}</div></td>
			<td class="entropy">${data.bwtsMtfEntropy}</td>
		</tr>`;
	}
	
	
	template += `
</table>`;

	resultsDiv.innerHTML = template;
}