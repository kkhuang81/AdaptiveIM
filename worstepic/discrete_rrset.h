/*
* BFS starting from one node
*/

int BuildHypergraphNode(int hyperId)
{
	const auto Idx = sfmt_genrand_uint32(&sfmtSeed) % NumcurNode;
	auto uStart = curNode[Idx];

	unsigned int n_visit_mark = 0, curIdx = 0;
	visit_mark[n_visit_mark++] = uStart;
	visit[uStart] = true;
	hyperG[uStart].push_back(hyperId);

	if (influModel == IC)
	{
		while (curIdx < n_visit_mark)
		{
			int i = visit_mark[curIdx++];
			for (int j = 0; j < (int)gT[i].size(); j++)
			{
				int v = gT[i][j];
				if (active_set[v] || visit[v])continue;
				double randDouble = sfmt_genrand_real1(&sfmtSeed);
				if (randDouble > probT[i][j])continue;

				visit[v] = true;
				visit_mark[n_visit_mark++] = v;
				hyperG[v].push_back(hyperId);
			}
		}
	}
	else if (influModel == LT)
	{
		while (curIdx < n_visit_mark)
		{
			int expand = visit_mark[curIdx++];

			if (gT[expand].size() == 0)
				continue;

			int index = sfmt_genrand_uint32(&sfmtSeed) % gT[expand].size();
			int v = gT[expand][index];
			if (active_set[v] || visit[v])continue;

			visit[v] = true;
			visit_mark[n_visit_mark++] = v;
			hyperG[v].push_back(hyperId);
		}
	}
	else
		ASSERT(false);

	for (unsigned int i = 0; i < n_visit_mark; i++)visit[visit_mark[i]] = false;

	hyperGT.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + n_visit_mark));

	return 0;
}

int BuildHypergraph(int hyperId)
{
	const auto Idx = sfmt_genrand_uint32(&sfmtSeed) % NumcurNode;
	auto uStart = curNode[Idx];

	unsigned int n_visit_mark = 0, curIdx = 0;
	visit_mark[n_visit_mark++] = uStart;
	visit[uStart] = true;
	hyperG_2[uStart].push_back(hyperId);

	if (influModel == IC)
	{
		while (curIdx < n_visit_mark)
		{
			int i = visit_mark[curIdx++];
			for (int j = 0; j < (int)gT[i].size(); j++)
			{
				int v = gT[i][j];
				if (active_set[v] || visit[v])continue;

				double randDouble = sfmt_genrand_real1(&sfmtSeed);
				if (randDouble > probT[i][j])continue;

				visit[v] = true;
				visit_mark[n_visit_mark++] = v;
				hyperG_2[v].push_back(hyperId);
			}
		}
	}
	else if (influModel == LT)
	{
		while (curIdx < n_visit_mark)
		{
			int expand = visit_mark[curIdx++];

			if (gT[expand].size() == 0)
				continue;

			int index = sfmt_genrand_uint32(&sfmtSeed) % gT[expand].size();
			int v = gT[expand][index];
			if (active_set[v] || visit[v])continue;

			visit[v] = true;
			visit_mark[n_visit_mark++] = v;
			hyperG_2[v].push_back(hyperId);
		}
	}
	else
		ASSERT(false);

	for (unsigned int i = 0; i < n_visit_mark; i++)visit[visit_mark[i]] = false;

	hyperGT_2.push_back(vector<int>(visit_mark.begin(), visit_mark.begin() + n_visit_mark));

	return 0;
}

//this function takes a quite long time
void realization(vector<int>& batch_set)
{
	unsigned int n_visit_mark = 0, curIdx = 0;
	for (auto seed : batch_set)
	{
		active_set[seed] = true;
		visit_mark[n_visit_mark++] = seed;
		--NumcurNode;
		curNodeIdx[curNode[NumcurNode]] = curNodeIdx[seed];
		curNode[curNodeIdx[seed]] = curNode[NumcurNode];
	}

	while (curIdx < n_visit_mark)
	{
		int expand = visit_mark[curIdx++];
		for (auto v : PO[expand])
		{
			if (active_set[v])continue;
			visit_mark[n_visit_mark++] = v;
			active_set[v] = true;
			--NumcurNode;
			curNodeIdx[curNode[NumcurNode]] = curNodeIdx[v];
			curNode[curNodeIdx[v]] = curNode[NumcurNode];
		}
	}

	//remove all RR sets

	for (auto& hyper : hyperG)hyper.clear();
	for (auto& hyperT : hyperGT)vector<int>().swap(hyperT);
	hyperGT.clear();

	for (auto& hyper : hyperG_2)hyper.clear();
	for (auto& hyperT : hyperGT_2)vector<int>().swap(hyperT);
	hyperGT_2.clear();

}
