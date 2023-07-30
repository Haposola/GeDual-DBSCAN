#pragma once
template<typename eleType>
class scoreStack {
public:
	scoreStack() {
		top = 0;
		//scoreArray.reserve(10);
		elementArray.reserve(10);
		//the initial size 10 is hard-coded. 
		//For dataset of size around 1000000, it usually need not to reallocate the storage of the stack.
	}

	inline void push( const eleType element) {
		/*if (top >= scoreArray.capacity()) {
			
			
			scoreArray.reserve(2 * scoreArray.capacity());
			elementArray.reserve(2 * scoreArray.capacity());
		}
		scoreArray[top] = score;
		elementArray[top] = element;*/
		//scoreArray.push_back(score);
		elementArray.push_back(element);
		top++;
	}
	inline void pop() {
		//scoreArray.pop_back();
		elementArray.pop_back();
		top--;
	}
	/*inline double topScore() {
		assert(top > 0);
		return scoreArray[top-1];
	}*/
	inline eleType topElement() {
		assert(top > 0);
		return elementArray[top - 1];
	}

	inline bool isEmpty() {
		return top == 0;
	}
	inline void clear() {
		//scoreArray.resize(0);
		elementArray.resize(0);
		//scoreArray.reserve(10);
		elementArray.reserve(10);
		top = 0;
	}


private:
	size_t top;
	//std::vector<double> scoreArray;
	std::vector<eleType> elementArray;
};