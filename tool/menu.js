//Struct MenuProperties
function MenuProperties() {
	//Options
	this.entropic = true;	
	this.nearEntropic = true;	
	this.bwts = true;	
	this.kSetShaping = 0;    
}
//End struct MenuProperties

//Class MenuManager
export function MenuManager() {
	
	this.properties = new MenuProperties();
	this.rootMenu = new dat.GUI();	
	
	//Options
	var optionsMenu = this.rootMenu.addFolder('Options');			
	optionsMenu.add(this.properties, 'entropic');		
	optionsMenu.add(this.properties, 'nearEntropic').name('near entropic');		
	optionsMenu.add(this.properties, 'bwts');								
	optionsMenu.add(this.properties, 'kSetShaping').name('k set shaping');				
	
}


var menuManager = new MenuManager();
window.menu = menuManager.properties;