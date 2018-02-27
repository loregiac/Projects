package jdraw.figures;

import java.awt.Point;
import java.awt.event.MouseEvent;

import jdraw.framework.DrawContext;
import jdraw.framework.Figure;
import jdraw.std.AddFigureCommand;

public abstract class AbstractDragDrawTool extends AbstractDrawTool {

	private Figure figure;
	private Point anchor;
	//private final DrawContext context = null;
	private DrawContext context;

	protected AbstractDragDrawTool(DrawContext context, String name, String icon) {
		super(name, icon);
		this.context = context;
	}
	
	
	
	/**
	 * Deactivates the current mode by resetting the cursor
	 * and clearing the status bar.
	 * @see jdraw.framework.DrawTool#deactivate()
	 */
	@Override
	public void deactivate() {
		context.showStatusText("");
	}

	/**
	 * Activates the Rectangle Mode. There will be a
	 * specific menu added to the menu bar that provides settings for
	 * Rectangle attributes
	 */
	@Override
	public void activate() {
		context.showStatusText(getName() + "Mode");
	}
	
	protected abstract Figure createFigure(Point p);
	
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
	public void mouseDown(int x, int y, MouseEvent e) {
		if (figure != null) {
			throw new IllegalStateException();
		}
		anchor = new Point(x, y);
		figure = createFigure(anchor);
		context.getModel().addFigure(figure);
	}

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
	public void mouseDrag(int x, int y, MouseEvent e) {
		figure.setBounds(anchor, new Point(x, y));
		java.awt.Rectangle r = figure.getBounds();
		this.context.showStatusText("w: " + r.width + ", h: " + r.height);
	}
	/**
	 * When the user releases the mouse, the Rectangle object is updated
	 * according to the color and fill status settings.
	 * 
	 * @param x   x-coordinate of mouse
	 * @param y   y-coordinate of mouse
	 * @param e   event containing additional information about which keys were
	 *            pressed.
	 * 
	 * @see jdraw.framework.DrawTool#mouseUp(int, int, MouseEvent)
	 */
	@Override
	public void mouseUp(int x, int y, MouseEvent e) {
		java.awt.Rectangle r = figure.getBounds();
		context.getModel().getDrawCommandHandler().addCommand(new AddFigureCommand(figure,context.getModel()));
		if (r.width == 0 && r.height == 0) {
			  context.getModel().removeFigure(figure); 
		}
		figure = null;
		anchor = null;
	}

	//@Override
	//public Cursor getCursor() {
		//return Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR); 
	//}
	
	

}
