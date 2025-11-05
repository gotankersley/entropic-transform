window.showResults = function(data) {	
	var resultsDiv = document.getElementById('results');
	var menu = window.menu;
	var template = 
`<table border="1">
	<tr>
		<th>Type</th>
		<th>Sequence</th>
		<th>Bijective Rank</th>
		<th>Entropy</th>
	</tr>
	<tr>
		<td>Random</td>
		<td><div class="seq">${data.rndSeq}</div></td>
		<td><div class="rank">${data.rndRank}</div></td>
		<td class="entropy">${data.rndEntropy}</td>
	</tr>`;
	
	if (menu.entropic) {
		template += 
		`<tr>
			<td>Entropic</td>
			<td><div class="seq">${data.entSeq}</div></td>
			<td><div class="rank">${data.rndRank}</div></td>
			<td class="entropy">${data.entEntropy}</td>
		</tr>`;
	
	}
	
	if (menu.nearEntropic) {
		template += 
		`<tr>
			<td>Near Entropic</td>
			<td><div class="seq">${data.nearSeq}</div></td>
			<td><div class="rank">${data.rndRank}</div></td>
			<td class="entropy">${data.nearEntropy}</td>
		</tr>`;
	}
	
	if (menu.bwts) {
		template += 
		`<tr>
			<td>BWTS</td>
			<td><div class="seq">${data.bwtsSeq}</div></td>
			<td><div class="rank">${data.bwtsRank}</div></td>
			<td class="entropy">${data.bwtsEntropy}</td>
		</tr>
		
		
		<tr>
			<td>BWTS + MTF </td>
			<td><div class="seq">${data.bwtsMtfSeq}</div></td>
			<td><div class="rank">${data.bwtsMtfRank}</div></td>
			<td class="entropy">${data.bwtsMtfEntropy}</td>
		</tr>`;
	}
	
	
	template += `
</table>`;

	resultsDiv.innerHTML = template;
}