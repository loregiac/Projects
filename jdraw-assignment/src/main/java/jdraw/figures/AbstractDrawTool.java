package jdraw.figures;

import java.awt.Cursor;
import java.awt.event.MouseEvent;

import javax.swing.Icon;
import javax.swing.ImageIcon;

import jdraw.framework.DrawTool;

public abstract class AbstractDrawTool implements DrawTool{
	private static final String IMAGES= "/images/";
	private String name;
	private String icon;
	
	protected AbstractDrawTool(String name, String icon){
		this.name = name;
		this.icon = icon;
	}
	
	@Override
	public Icon getIcon() {
		if(icon!=null){
			return new ImageIcon(getClass().getResource(IMAGES + icon));
		}else{
			return null;
		}
	}

	@Override
	public final String getName() {
		return name;
	}
	
	/**
	 * Deactivates the current mode by resetting the cursor
	 * and clearing the status bar.
	 * @see jdraw.framework.DrawTool#deactivate()
	 */
	@Override
	public void deactivate() {}

	/**
	 * Activates the Rectangle Mode. There will be a
	 * specific menu added to the menu bar that provides settings for
	 * Rectangle attributes
	 */
	@Override
	public void activate(){};
	
	/**
	 * Initializes a new Rectangle object by setting an anchor
	 * point where the mouse was pressed. A new Rectangle is then
	 * added to the model.
	 * @param x x-coordinate of mouse
	 * @param y y-coordinate of mouse
	 * @param e event containing additional information about which keys were pressed.
	 * 
	 * @see jdraw.framework.DrawTool#mouseDown(int, int, MouseEvent)
	 */
	@Override
	public void mouseDown(int x, int y, MouseEvent e) {}

	/**
	 * During a mouse drag, the Rectangle will be resized according to the mouse
	 * position. The status bar shows the current size.
	 * 
	 * @param x   x-coordinate of mouse
	 * @param y   y-coordinate of mouse
	 * @param e   event containing additional information about which keys were
	 *            pressed.
	 * 
	 * @see jdraw.framework.DrawTool#mouseDrag(int, int, MouseEvent)
	 */
	@Override
	public void mouseDrag(int x, int y, MouseEvent e) {}
	
	@Override
	public void mouseUp(int x, int y, MouseEvent e) {}
	
	@Override
	   public Cursor getCursor() {
	return Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR); }

	
}
