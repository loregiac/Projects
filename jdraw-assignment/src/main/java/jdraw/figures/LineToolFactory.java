package jdraw.figures;

import jdraw.framework.DrawContext;
import jdraw.framework.DrawTool;
import jdraw.framework.DrawToolFactory;

public class LineToolFactory implements DrawToolFactory {

	private String name;
	
	private String iconName;

	
	@Override
	public String getName() {
		
		return name;
	}

	@Override
	public void setName(String name) {
		this.name = name;
	}

	@Override
	public String getIconName() {
		// TODO Auto-generated method stub
		return iconName;
	}

	@Override
	public void setIconName(String name) {
		// TODO Auto-generated method stub
		this.iconName=name;
	}

	@Override
	public DrawTool createTool(DrawContext context) {
		return new LineTool(context,name,iconName);
	}

}
