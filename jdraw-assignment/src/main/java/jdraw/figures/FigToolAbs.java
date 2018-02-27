package jdraw.figures;

import java.awt.Cursor;

import jdraw.framework.DrawContext;
import jdraw.framework.DrawView;

public abstract class FigToolAbs extends FigTool{

	/**
	 * The context we use for drawing.
	 */
	protected DrawContext context;
	/**
	 * The context's view. This variable can be used as a shortcut, i.e.
	 * instead of calling context.getView().
	 */
	protected DrawView view;
	/**
	 * The context's view. This variable can be used as a shortcut, i.e.
	 * instead of calling context.getView().
	 */
	
	public FigToolAbs(DrawContext context) {
		super(context);
	}
	
	@Override
	public void deactivate() {
		this.context.showStatusText("");
	}

	@Override
	public Cursor getCursor() {
		return Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR);
	}
	
}
