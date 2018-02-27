/*
 * Copyright (c) 2017 Fachhochschule Nordwestschweiz (FHNW)
 * All Rights Reserved.
 */
package jdraw.std;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.LinkedList;
import java.util.List;

import javax.swing.ButtonGroup;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;
import javax.swing.filechooser.FileFilter;

import jdraw.Grids.*;
import jdraw.figures.FixTool;
import jdraw.figures.LineTool;
import jdraw.figures.OvalTool;
import jdraw.figures.RectTool;
import jdraw.framework.DrawCommandHandler;
import jdraw.framework.DrawModel;
import jdraw.framework.DrawTool;
import jdraw.framework.DrawToolFactory;
import jdraw.framework.DrawView;
import jdraw.framework.Figure;
import jdraw.actions.*;
/**
 * Standard implementation of interface DrawContext.
 * 
 * @see DrawView
 * @author Dominik Gruntz & Christoph Denzler
 * @version 2.6, 24.09.09
 */
public class StdContext extends AbstractContext {

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	/**
	 * Constructs a standard context with a default set of drawing tools.
	 * @param view the view that is displaying the actual drawing.
	 */
  public StdContext(DrawView view) {
		super(view, null);
	}
	
  /**
   * Constructs a standard context. The drawing tools available can be parameterized using <code>toolFactories</code>.
   * @param view the view that is displaying the actual drawing.
   * @param toolFactories a list of DrawToolFactories that are available to the user
   */
	public StdContext(DrawView view, List<DrawToolFactory> toolFactories) {
		super(view, toolFactories);
	}

	/**
	 * Creates and initializes the "Edit" menu.
	 * 
	 * @return the new "Edit" menu.
	 */
	@Override
	protected JMenu createEditMenu() {
		JMenu editMenu = new JMenu("Edit");
		final JMenuItem undo = new JMenuItem("Undo");
		undo.setAccelerator(KeyStroke.getKeyStroke("control Z"));
		editMenu.add(undo);
		undo.addActionListener(e -> {
				final DrawCommandHandler h = getModel().getDrawCommandHandler();
				if (h.undoPossible()) {
					h.undo();
				}
			}
		);

		final JMenuItem redo = new JMenuItem("Redo");
		redo.setAccelerator(KeyStroke.getKeyStroke("control Y"));
		editMenu.add(redo);
		redo.addActionListener(e -> {
				final DrawCommandHandler h = getModel().getDrawCommandHandler();
				if (h.redoPossible()) {
					h.redo();
				}
			}
		);
		editMenu.addSeparator();

		JMenuItem sa = new JMenuItem("SelectAll");
		sa.setAccelerator(KeyStroke.getKeyStroke("control A"));
		editMenu.add(sa);
		sa.addActionListener( e -> {
				for (Figure f : getModel().getFigures()) {
					getView().addToSelection(f);
				}
				getView().repaint();
			}
		);

		editMenu.addSeparator();
		//editMenu.add("Cut").setEnabled(false);
		//editMenu.add(new CopyAction(this, editMenu));
		JMenuItem copy = new JMenuItem(new CopyAction(this, editMenu));
		copy.setAccelerator(KeyStroke.getKeyStroke("control C"));
		editMenu.add(copy);
		JMenuItem cut = new JMenuItem(new CutAction(this, editMenu));
		cut.setAccelerator(KeyStroke.getKeyStroke("control X"));
		editMenu.add(cut);
		
		//editMenu.add("Copy").setEnabled(false);
		//editMenu.add("Paste").setEnabled(false);
		//editMenu.add(new PasteAction(this, editMenu));
		JMenuItem paste = new JMenuItem(new PasteAction(this, editMenu));
		paste.setAccelerator(KeyStroke.getKeyStroke("control V"));
		editMenu.add(paste);
		
		editMenu.addSeparator();

		JMenuItem borderDecorator = new JMenuItem(new BorderDecoratorAction(this, editMenu));
		editMenu.add(borderDecorator);
		
		JMenuItem bundleDecorator = new JMenuItem(new BundleDecoratorAction(this, editMenu));
		editMenu.add(bundleDecorator);

		JMenuItem animationDecorator = new JMenuItem(new AnimationDecoratorAction(this, editMenu));
		editMenu.add(animationDecorator);
		
		editMenu.addSeparator();
		JMenuItem clear = new JMenuItem("Clear");
		editMenu.add(clear);
		clear.addActionListener(e -> {
			getModel().removeAllFigures();
		});
		
		editMenu.addSeparator();
		//JMenuItem group = new JMenuItem("Group");
		//group.setEnabled(false);
		//editMenu.add(group);
		editMenu.add(new GroupAction(this, editMenu));
		
		//JMenuItem ungroup = new JMenuItem("Ungroup");
		//ungroup.setEnabled(false);
		//editMenu.add(ungroup);
		editMenu.add(new UngroupAction(this, editMenu));

		
		editMenu.addSeparator();

		JMenu orderMenu = new JMenu("Order...");
		JMenuItem frontItem = new JMenuItem("Bring To Front");
		frontItem.addActionListener(e -> {
			bringToFront(getView().getModel(), getView().getSelection());
		});
		orderMenu.add(frontItem);
		JMenuItem backItem = new JMenuItem("Send To Back");
		backItem.addActionListener(e -> {
			sendToBack(getView().getModel(), getView().getSelection());
		});
		orderMenu.add(backItem);
		editMenu.add(orderMenu);
		
		editMenu.addSeparator();

		JMenuItem fix = new JMenuItem("Fix");
		editMenu.add(fix);
		FixTool fixTool = new FixTool();
		fix.addActionListener(new ActionListener(){
			@Override
			public void actionPerformed(ActionEvent e) {
				List<Figure> list = getView().getSelection();
				for(Figure f : list){
					fixTool.toggleFixFigure(f);
				}
			}
		});
		
		JMenu grid = new JMenu("Grid...");
		
		JCheckBoxMenuItem noGrid = new JCheckBoxMenuItem("No Grid");
		noGrid.addActionListener(e -> {
			super.getView().setConstrainer(null);
		});
		
		JCheckBoxMenuItem grid10 = new JCheckBoxMenuItem("Grid 10");
		grid10.addActionListener(e -> {
			super.getView().setConstrainer(new StepGrid(10,10));
		});
		
		JCheckBoxMenuItem grid20 = new JCheckBoxMenuItem("Grid 20");
		grid20.addActionListener(e -> {
			super.getView().setConstrainer(new StepGrid(20,20));
		});
		
		JCheckBoxMenuItem grid50 = new JCheckBoxMenuItem("Grid 50");
		grid50.addActionListener(e -> {
			super.getView().setConstrainer(new StepGrid(50,50));
		});
		
	    ButtonGroup bg = new ButtonGroup();
	    bg.add(noGrid);
	    bg.add(grid10);
		bg.add(grid20);
		bg.add(grid50);
		noGrid.setSelected(true);
		
	    grid.add(noGrid);
		grid.add(grid10);
		grid.add(grid20);
		grid.add(grid50);
		editMenu.add(grid);
		
		return editMenu;
	}



	/**
	 * Creates and initializes items in the file menu.
	 * 
	 * @return the new "File" menu.
	 */
	@Override
	protected JMenu createFileMenu() {
	  JMenu fileMenu = new JMenu("File");
		JMenuItem open = new JMenuItem("Open");
		fileMenu.add(open);
		open.setAccelerator(KeyStroke.getKeyStroke("control O"));
		open.addActionListener(e -> doOpen());

		JMenuItem save = new JMenuItem("Save");
		save.setAccelerator(KeyStroke.getKeyStroke("control S"));
		fileMenu.add(save);
		save.addActionListener(e ->	doSave());

		JMenuItem exit = new JMenuItem("Exit");
		fileMenu.add(exit);
		exit.addActionListener(e -> System.exit(0));
		
		return fileMenu;
	}

	@Override
	protected void doRegisterDrawTools() {
		
		/*DrawTool rectangleTool = new RectTool(this,"Rectangle","rectangle.png");
		DrawTool ovalTool = new OvalTool(this,"Oval","oval.png");
		DrawTool lineTool = new LineTool(this,"Line","line.png");

		addTool(rectangleTool);
		addTool(ovalTool);
		addTool(lineTool);*/
		for(DrawToolFactory d : super.getToolFactories()){
			addTool(d.createTool(this));
		}

	}

	/**
	 * Changes the order of figures and moves the figures in the selection
	 * to the front, i.e. moves them to the end of the list of figures.
	 * @param model model in which the order has to be changed
	 * @param selection selection which is moved to front
	 */
	public void bringToFront(DrawModel model, List<Figure> selection) {
		// the figures in the selection are ordered according to the order in
		// the model
		List<Figure> orderedSelection = new LinkedList<Figure>();
		int pos = 0;
		for (Figure f : model.getFigures()) {
			pos++;
			if (selection.contains(f)) {
				orderedSelection.add(0, f);
			}
		}
		for (Figure f : orderedSelection) {
			model.setFigureIndex(f, --pos);
		}
	}

	/**
	 * Changes the order of figures and moves the figures in the selection
	 * to the back, i.e. moves them to the front of the list of figures.
	 * @param model model in which the order has to be changed
	 * @param selection selection which is moved to the back
	 */
	public void sendToBack(DrawModel model, List<Figure> selection) {
		// the figures in the selection are ordered according to the order in
		// the model
		List<Figure> orderedSelection = new LinkedList<Figure>();
		for (Figure f : model.getFigures()) {
			if (selection.contains(f)) {
				orderedSelection.add(f);
			}
		}
		int pos = 0;
		for (Figure f : orderedSelection) {
			model.setFigureIndex(f, pos++);
		}
	}

	/**
	 * Handles the saving of a drawing to a file.
	 */
	private void doSave() {
		JFileChooser chooser = new JFileChooser(getClass().getResource("")
				.getFile());
		chooser.setDialogTitle("Save Graphic");
		chooser.setDialogType(JFileChooser.SAVE_DIALOG);
		FileFilter filter = new FileFilter() {
			@Override
			public String getDescription() {
				return "JDraw Graphic (*.draw)";
			}

			@Override
			public boolean accept(File f) {
				return f.getName().endsWith(".draw");
			}
		};
		chooser.setFileFilter(filter);
		int res = chooser.showSaveDialog(this);

		if (res == JFileChooser.APPROVE_OPTION) {
			// save graphic
			File file = chooser.getSelectedFile();
			if (chooser.getFileFilter() == filter && !filter.accept(file)) {
				file = new File(chooser.getCurrentDirectory(), file.getName() + ".draw");
			}
			System.out.println("save current graphic to file " + file.getName());
		}
	}

	/**
	 * Handles the opening of a new drawing from a file.
	 */
	private void doOpen() {
		JFileChooser chooser = new JFileChooser(getClass().getResource("")
				.getFile());
		chooser.setDialogTitle("Open Graphic");
		chooser.setDialogType(JFileChooser.OPEN_DIALOG);
		chooser.setFileFilter(new javax.swing.filechooser.FileFilter() {
			@Override
			public String getDescription() {
				return "JDraw Graphic (*.draw)";
			}

			@Override
			public boolean accept(File f) {
				return f.isDirectory() || f.getName().endsWith(".draw");
			}
		});
		int res = chooser.showOpenDialog(this);

		if (res == JFileChooser.APPROVE_OPTION) {
			// read jdraw graphic
			System.out.println("read file "
					+ chooser.getSelectedFile().getName());
		}
	}

}
