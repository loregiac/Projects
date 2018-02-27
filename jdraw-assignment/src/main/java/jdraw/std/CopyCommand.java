/*
 * Copyright (c) 2017 Fachhochschule Nordwestschweiz (FHNW)
 * All Rights Reserved. 
 */

package jdraw.std;

import jdraw.framework.DrawCommand;
import jdraw.framework.Figure;

public class CopyCommand implements DrawCommand {
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;


	private Figure clone;

	public CopyCommand(Figure clone) {

		this.clone = clone;
	}

	@Override
	public void redo() {
		ClipBoard.add(clone);
	}

	@Override
	public void undo() {
		ClipBoard.remove(clone);
	}

}

